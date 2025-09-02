"""

Pełne GUI - łączy moduł generator/generacja.py z interfejsem Tkinter,
rysuje wykres czasowy i FFT, obsługuje manual_sweep (w tle), synchronizuje
dostęp do oscyloskopu lockiem (osc_lock).
"""

import os
import sys
import threading
import time
import math
import tkinter as tk
from tkinter import ttk, messagebox

import numpy as np
import pyvisa
import matplotlib
# backend 'TkAgg' potrzebny do interaktywnego GUI; jeśli masz headless server użyj 'Agg'
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# dodaj ścieżkę do folderu generator jeżeli katalog inny
script_dir = os.path.dirname(os.path.abspath(__file__))
gen_dir = os.path.join(script_dir, "generator")
if gen_dir not in sys.path:
    sys.path.append(gen_dir)

from generacja import init_device  # import klasy z generacja.py

# --------------------- KONFIG: ustaw tu adresy swoich urządzeń ---------------------
GEN_ADDR = "TCPIP::192.168.0.200::INSTR"   # generator
OSC_ADDR = "TCPIP0::192.168.0.231::inst0::INSTR"  # oscyloskop (dostosuj: format może zależeć od sterownika)
# -------------------------------------------------------------------------------

# czas pomiędzy odświeżeniami wykresu (s)
PLOT_INTERVAL = 0.5


# --------------------- pomocnicze funkcje ---------------------
def to_seconds(val_str, unit):
    if val_str is None or val_str == "":
        return None
    raw = float(val_str)
    mul = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "ns": 1e-9}.get(unit, 1.0)
    return raw * mul

def to_hz(val_str, unit):
    if val_str is None or val_str == "":
        return None
    raw = float(val_str)
    mul = {"MHz": 1e6, "kHz": 1e3, "Hz": 1.0, "mHz": 1e-3, "µHz": 1e-6}.get(unit, 1.0)
    return raw * mul


# --------------------- APLIKACJA ---------------------
class WielkiBratApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Wielki Brat Premium")
        self.root.geometry("1100x820")

        # device: generator + oscilloscope
        self.device = init_device(generator_address=GEN_ADDR, oscilloscope_address=OSC_ADDR)
        self.osc_lock = threading.Lock()  # lock do synchronizacji VISA
        # (przekażemy ten lock do manual_sweep aby update_plot nie walczył z nim)

        # spróbuj połączyć (nie powoduj wyjątku - obsłuż komunikatem)
        try:
            idn_gen = self.device.connect_generator(GEN_ADDR)
            print("Generator IDN:", idn_gen)
        except Exception as e:
            print("Generator connect warning:", e)
        try:
            idn_osc = self.device.connect_oscilloscope(OSC_ADDR)
            print("Oscilloscope IDN:", idn_osc)
        except Exception as e:
            print("Oscilloscope connect warning:", e)

        # domyślne skalowania odczytane z oscyloskopu (jeśli podłączony)
        try:
            if self.device.oscilloscope:
                self.ymult = float(self.device.oscilloscope.query("WFMPRE:YMULT?"))
                self.yzero = float(self.device.oscilloscope.query("WFMPRE:YZERO?"))
                self.yoff  = float(self.device.oscilloscope.query("WFMPRE:YOFF?"))
                self.xincr = float(self.device.oscilloscope.query("WFMPRE:XINCR?"))
            else:
                # bez oscyloskopu użyj sensownych wartości
                self.ymult = 1.0
                self.yzero = 0.0
                self.yoff = 0.0
                self.xincr = 1e-6
        except Exception as e:
            print("Nie udało się pobrać WFMPRE (używam domyślnych):", e)
            self.ymult = 1.0
            self.yzero = 0.0
            self.yoff = 0.0
            self.xincr = 1e-6

        # zmienne GUI
        self.channel_var = tk.StringVar(value="C1")
        self.osc_channel_var = tk.StringVar(value="CH1")
        self.vs_var = tk.StringVar()   # V/div
        self.hs_var = tk.StringVar()   # czas/div - liczba
        self.hs_unit_var = tk.StringVar(value="s")  # jednostka czasu
        self.current_load = tk.StringVar(value="50 Ω")

        # budowa GUI
        self.build_gui()

        # tworzymy figure 2 subplots (time + FFT)
        self.fig, (self.ax_time, self.ax_freq) = plt.subplots(2, 1, figsize=(7, 8))
        self.fig.tight_layout(pad=3.0)
        self.line_time, = self.ax_time.plot([], [], lw=1.5)
        self.line_fft,  = self.ax_freq.plot([], [], lw=1.0)

        self.ax_time.set_title("Sygnał w czasie rzeczywistym")
        self.ax_time.set_xlabel("Czas [s]")
        self.ax_time.set_ylabel("Napięcie [V]")
        self.ax_time.grid(True)

        self.ax_freq.set_title("Widmo amplitudowe (FFT)")
        self.ax_freq.set_xlabel("Częstotliwość [Hz]")
        self.ax_freq.set_ylabel("Amplituda")
        self.ax_freq.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=8, pady=8)

        # wątek aktualizacji wykresu
        self._stop_plot_thread = False
        self.plot_thread = threading.Thread(target=self.update_plot_loop, daemon=True)
        self.plot_thread.start()

    def build_gui(self):
        # layout: left control panel + right plot panel
        container = tk.Frame(self.root)
        container.pack(fill="both", expand=True)

        left = tk.Frame(container, width=420, bg="lightgray", padx=10, pady=10)
        left.pack(side="left", fill="y")
        self.right_frame = tk.Frame(container, bg="white")
        self.right_frame.pack(side="right", fill="both", expand=True)

        # --------- Left panel controls ----------
        ttk.Label(left, text="Kanał generatora:", background="lightgray").grid(row=0, column=0, sticky="w")
        ttk.Combobox(left, textvariable=self.channel_var, values=["C1","C2"], state="readonly", width=8).grid(row=0, column=1, sticky="w")

        ttk.Label(left, text="Kanał oscyloskopu:", background="lightgray").grid(row=1, column=0, sticky="w", pady=(6,0))
        ttk.Combobox(left, textvariable=self.osc_channel_var,
                     values=["CH1","CH2","CH3","CH4"], state="readonly", width=8).grid(row=1, column=1, sticky="w", pady=(6,0))

        ttk.Label(left, text="Skala pionowa (V/div):", background="lightgray").grid(row=2, column=0, sticky="w", pady=(8,0))
        tk.Entry(left, textvariable=self.vs_var, justify="center", width=12).grid(row=2, column=1, sticky="w", pady=(8,0))

        ttk.Label(left, text="Skala pozioma (czas/div):", background="lightgray").grid(row=3, column=0, sticky="w", pady=(6,0))
        tk.Entry(left, textvariable=self.hs_var, justify="center", width=8).grid(row=3, column=1, sticky="w", pady=(6,0))
        ttk.Combobox(left, textvariable=self.hs_unit_var, values=["s","ms","us","ns"], state="readonly", width=6).grid(row=3, column=2, sticky="w", padx=(5,0), pady=(6,0))

        # LOAD
        ttk.Label(left, text="Ustaw impedancję LOAD:", background="lightgray").grid(row=4, column=0, sticky="w", pady=(8,0))
        self.load_mode_var = tk.StringVar(value="50")
        ttk.Combobox(left, textvariable=self.load_mode_var, values=["50","HZ","Własna"], state="readonly", width=8).grid(row=4, column=1, sticky="w")
        self.custom_load_var = tk.StringVar()
        self.custom_load_entry = tk.Entry(left, textvariable=self.custom_load_var, justify="center", width=8)

        def on_load_change(*args):
            mode = self.load_mode_var.get()
            ch = self.channel_var.get()
            if mode == "Własna":
                self.custom_load_entry.grid(row=4, column=2, sticky="w", padx=(5,0))
            else:
                self.custom_load_entry.grid_forget()
            # spróbuj wysłać do generatora
            try:
                if mode == "Własna":
                    val = float(self.custom_load_var.get())
                    self.device.SET_LOAD(channel=ch, load=val)
                    self.current_load.set(f"{val:.0f} Ω")
                else:
                    self.device.SET_LOAD(channel=ch, load=mode)
                    self.current_load.set("High Z" if mode == "HZ" else f"{mode} Ω")
            except Exception:
                pass

        self.load_mode_var.trace_add("write", on_load_change)
        self.custom_load_var.trace_add("write", on_load_change)

        # Wave combine
        self.wave_combine_var = tk.BooleanVar(value=False)
        tk.Checkbutton(left, text="✔ Wave Combine", bg="lightgray",
                       variable=self.wave_combine_var,
                       command=lambda: self.device.WAVE_COMBINE(channel=self.channel_var.get(), state=self.wave_combine_var.get()),
                       ).grid(row=5, column=0, columnspan=2, sticky="w", pady=(10,6))

        # signal type combobox + dynamic form area
        ttk.Label(left, text="Typ sygnału:", background="lightgray").grid(row=6, column=0, sticky="w", pady=(10,0))
        self.function_choices = ["SINE","SQUARE","RAMP","PULSE","NOISE","ARB","DC","MANUAL SWEEP"]
        self.func_cb = ttk.Combobox(left, values=self.function_choices, state="readonly", width=12)
        self.func_cb.grid(row=6, column=1, sticky="w", pady=(10,0))
        self.func_cb.set(self.function_choices[0])

        # dynamiczny formularz
        self.left_form = tk.Frame(left, bg="lightgray")
        self.left_form.grid(row=7, column=0, columnspan=3, sticky="nsew", pady=(10,0))

        # binder
        self.func_cb.bind("<<ComboboxSelected>>", lambda e: self.update_form(self.func_cb.get()))
        self.update_form(self.function_choices[0])

    # ---------- helper do tworzenia entry + unit ----------
    def create_entry(self, parent, label_text, width=25):
        wrapper = tk.Frame(parent, bg="lightgray")
        wrapper.pack(pady=5, fill="x")
        tk.Label(wrapper, text=label_text, bg="lightgray").pack()
        ent = tk.Entry(wrapper, justify="center", width=width)
        ent.pack(pady=(2,0))
        return ent

    def create_freq_entry_with_unit(self, parent, label_text, entries_dict, key_base, default_unit="Hz"):
        frame = tk.Frame(parent, bg="lightgray")
        frame.pack(pady=5, fill="x")
        tk.Label(frame, text=label_text, bg="lightgray").pack()
        inner = tk.Frame(frame, bg="lightgray"); inner.pack(pady=(2,0))
        var_val = tk.StringVar()
        var_unit = tk.StringVar(value=default_unit)
        tk.Entry(inner, textvariable=var_val, justify="center", width=12).pack(side="left")
        combo = ttk.Combobox(inner, textvariable=var_unit, values=["MHz","kHz","Hz","mHz","µHz"], state="readonly", width=5)
        combo.pack(side="left", padx=(5,0)); combo.set(default_unit)
        entries_dict[key_base + "_val"] = var_val
        entries_dict[key_base + "_unit"] = var_unit
        return var_val, var_unit

    # ---------- update_form (dynamiczny formularz) ----------
    def update_form(self, selected_function):
        # wyczyść poprzednie widgety
        for w in self.left_form.winfo_children():
            w.destroy()

        def create_entry_local(label):
            return self.create_entry(self.left_form, label)

        def create_freq_local(label, key):
            return self.create_freq_entry_with_unit(self.left_form, label, entries, key)

        entries = {}

        if selected_function in ["SINE","SQUARE","RAMP","ARB"]:
            # frequency with units
            create_freq_local("Frequency", "freq")
            entries["amp"] = create_entry_local("Amplitude [Vpp]")
            entries["offset"] = create_entry_local("Offset [Vdc]")
            entries["phase"] = create_entry_local("Phase")
            if selected_function == "SQUARE":
                entries["duty"] = create_entry_local("Duty [%]")
            if selected_function == "RAMP":
                entries["sym"] = create_entry_local("Symmetry [%]")

            # sweep controls
            self._create_sweep_options(entries)

            def call_fn():
                try:
                    freq_raw = entries["freq_val"].get().strip()
                    freq_unit = entries["freq_unit"].get()
                    frequency = to_hz(freq_raw, freq_unit) if freq_raw else 0.0
                    amp = float(entries["amp"].get())
                    off = float(entries["offset"].get())
                    phase = float(entries["phase"].get() or 0)
                    ch = self.channel_var.get()
                    if selected_function == "SINE":
                        self.device.SINE(channel=ch, frequency=frequency, amplitude=amp, offset=off, phase=phase)
                    elif selected_function == "SQUARE":
                        duty = float(entries["duty"].get() or 50)
                        self.device.SQUARE(channel=ch, frequency=frequency, amplitude=amp, offset=off, phase=phase, duty=duty)
                    elif selected_function == "RAMP":
                        sym = float(entries["sym"].get() or 50)
                        self.device.RAMP(channel=ch, frequency=frequency, amplitude=amp, offset=off, phase=phase, symmetry=sym)
                    elif selected_function == "ARB":
                        self.device.ARB(channel=ch, frequency=frequency, amplitude=amp, offset=off, phase=phase)

                    # if sweep enabled -> call generator.SWEEP_hw (optional)
                    if hasattr(self, "sweep_enabled") and self.sweep_enabled.get():
                        try:
                            start = float(self.sweep_fields["start_freq"].get())
                            stop = float(self.sweep_fields["stop_freq"].get())
                            sweep_time = float(self.sweep_fields["sweep_time"].get())
                            mode = self.sweep_fields["mode"].get()
                            direction = self.sweep_fields["direction"].get()
                            self.device.SWEEP_hw(channel=ch, start_freq=start, stop_freq=stop,
                                                 sweep_time=sweep_time, mode=mode, direction=direction,
                                                 wave_type=selected_function, amplitude=amp, offset=off, frequency=frequency)
                        except Exception as e:
                            print("Sweep_hw config error:", e)
                except Exception as e:
                    messagebox.showerror("Błąd", str(e))

            tk.Button(self.left_form, text="Zatwierdź", command=call_fn).pack(pady=8)

        elif selected_function == "PULSE":
            create_freq_local("Frequency", "freq")
            entries["duty"] = create_entry_local("Duty [%]")
            entries["amp"] = create_entry_local("Amplitude")
            entries["offset"] = create_entry_local("Offset")
            entries["pulse_width"] = create_entry_local("Pulse Width [s]")
            entries["rise"] = create_entry_local("Rise Edge [s]")
            entries["fall"] = create_entry_local("Fall Edge [s]")
            entries["delay"] = create_entry_local("Delay [s]")

            def run_pulse():
                try:
                    raw = entries["freq_val"].get().strip()
                    unit = entries["freq_unit"].get()
                    frequency = to_hz(raw, unit) if raw else 0.0
                    self.device.PULSE(channel=self.channel_var.get(),
                                      frequency=frequency,
                                      duty=float(entries["duty"].get() or 50.0),
                                      amplitude=float(entries["amp"].get()),
                                      offset=float(entries["offset"].get()),
                                      pulse_width=float(entries["pulse_width"].get()),
                                      rise_edge=float(entries["rise"].get()),
                                      fall=float(entries["fall"].get()),
                                      delay=float(entries["delay"].get()))
                except Exception as e:
                    messagebox.showerror("Błąd", str(e))

            tk.Button(self.left_form, text="Zatwierdź", command=run_pulse).pack(pady=8)

        elif selected_function == "NOISE":
            entries["stdev"] = create_entry_local("Std Dev")
            entries["mean"] = create_entry_local("Mean")
            def run_noise():
                try:
                    self.device.NOISE(channel=self.channel_var.get(),
                                      stdev=float(entries["stdev"].get()),
                                      mean=float(entries["mean"].get()))
                except Exception as e:
                    messagebox.showerror("Błąd", str(e))
            tk.Button(self.left_form, text="Zatwierdź", command=run_noise).pack(pady=8)

        elif selected_function == "DC":
            entries["offset"] = create_entry_local("Offset")
            def run_dc():
                try:
                    self.device.DC(channel=self.channel_var.get(), offset=float(entries["offset"].get()))
                except Exception as e:
                    messagebox.showerror("Błąd", str(e))
            tk.Button(self.left_form, text="Zatwierdź", command=run_dc).pack(pady=8)

        elif selected_function == "MANUAL SWEEP":
            # pola Start/Stop/Step z jednostkami
            create_freq_local("Start freq", "start")
            create_freq_local("Stop freq", "stop")
            create_freq_local("Step", "step")

            entries["amplitude"] = create_entry_local("Amplitude [Vpp]")
            entries["offset"] = create_entry_local("Offset [Vdc]")
            entries["delay"] = create_entry_local("Delay [s]")

            f = tk.Frame(self.left_form, bg="lightgray"); f.pack(pady=5, fill="x")
            tk.Label(f, text="Typ fali (SINE,SQUARE,RAMP,ARB,DC,NOISE,PULSE)", bg="lightgray").pack()
            wave_type_entry = tk.Entry(f, justify="center", width=25); wave_type_entry.insert(0, "SINE"); wave_type_entry.pack(pady=(2,0))

            # przycisk uruchamiający sweep w tle
            btn = tk.Button(self.left_form, text="Zatwierdź Sweep")
            btn.pack(pady=(10,20))

            def run_manual_sweep():
                # parse freq -> Hz
                def to_hz_local(raw, unit):
                    mapping = {"MHz":1e6, "kHz":1e3, "Hz":1.0, "mHz":1e-3, "µHz":1e-6}
                    return float(raw) * mapping[unit]

                try:
                    raw_start = entries["start_val"].get().strip()
                    raw_stop = entries["stop_val"].get().strip()
                    raw_step = entries["step_val"].get().strip()
                    if raw_start == "" or raw_stop == "" or raw_step == "":
                        messagebox.showerror("Błąd", "Uzupełnij pola Start/Stop/Step")
                        return

                    start_hz = to_hz_local(raw_start, entries["start_unit"].get())
                    stop_hz  = to_hz_local(raw_stop, entries["stop_unit"].get())
                    step_hz  = to_hz_local(raw_step, entries["step_unit"].get())
                    amp = float(entries["amplitude"].get())
                    off = float(entries["offset"].get())
                    dly = float(entries["delay"].get())
                    wav = wave_type_entry.get().strip().upper()
                except Exception as e:
                    messagebox.showerror("Błąd parsowania", str(e))
                    return

                # worker (wątek)
                def worker():
                    btn.config(state="disabled", text="Sweep trwa…")
                    try:
                        # użyj manual_sweep z urządzenia, przekaż lock by synchronizować I/O
                        self.device.manual_sweep(
                            channel=self.channel_var.get(),
                            start_freq=start_hz,
                            stop_freq=stop_hz,
                            step=step_hz,
                            amplitude=amp,
                            offset=off,
                            delay=dly,
                            oscilloscope=self.device.oscilloscope,
                            ymult=self.ymult, yzero=self.yzero, yoff=self.yoff, xincr=self.xincr,
                            osc_channel=self.osc_channel_var.get(),
                            csv_file="manual_sweep_output.csv",
                            wave_type=wav,
                            osc_lock_local=self.osc_lock,
                            progress_callback=None
                        )
                        messagebox.showinfo("Gotowe", "Manual sweep zakończony i zapisany do manual_sweep_output.csv")
                    except Exception as e:
                        messagebox.showerror("Błąd sweepu", str(e))
                    finally:
                        btn.config(state="normal", text="Zatwierdź Sweep")

                threading.Thread(target=worker, daemon=True).start()

            # attach
            btn.config(command=run_manual_sweep)

    # sweep helper UI
    def _create_sweep_options(self, entries_store):
        # pokazywanie prostego bloku sweep (używane w kilku miejscach)
        self.sweep_enabled = tk.BooleanVar(value=False)
        sweep_box = tk.Checkbutton(self.left_form, text="✔ Włącz sweep", variable=self.sweep_enabled, bg="lightgray")
        sweep_box.pack(pady=4)
        # fields dict (start/stop/time/mode/direction)
        self.sweep_fields = {}
        wrapper = tk.Frame(self.left_form, bg="lightgray"); wrapper.pack(pady=4, fill="x")
        tk.Label(wrapper, text="Sweep start [Hz]", bg="lightgray").pack()
        e1 = tk.Entry(wrapper, justify="center", width=25); e1.pack()
        self.sweep_fields["start_freq"] = e1
        tk.Label(wrapper, text="Sweep stop [Hz]", bg="lightgray").pack()
        e2 = tk.Entry(wrapper, justify="center", width=25); e2.pack()
        self.sweep_fields["stop_freq"] = e2
        tk.Label(wrapper, text="Sweep time [s]", bg="lightgray").pack()
        e3 = tk.Entry(wrapper, justify="center", width=25); e3.pack()
        self.sweep_fields["sweep_time"] = e3
        # mode/direction
        combo_frame = tk.Frame(self.left_form, bg="lightgray"); combo_frame.pack(pady=4)
        tk.Label(combo_frame, text="Sweep mode:", bg="lightgray").grid(row=0, column=0)
        cb1 = ttk.Combobox(combo_frame, values=["LINE","LOG"], state="readonly", width=10); cb1.grid(row=0, column=1)
        cb1.set("LINE")
        tk.Label(combo_frame, text="Direction:", bg="lightgray").grid(row=1, column=0)
        cb2 = ttk.Combobox(combo_frame, values=["UP","DOWN"], state="readonly", width=10); cb2.grid(row=1, column=1)
        cb2.set("UP")
        self.sweep_fields["mode"] = cb1; self.sweep_fields["direction"] = cb2

    # ---------- pętla aktualizacji wykresów ----------
    def update_plot_loop(self):
        while not self._stop_plot_thread:
            try:
                # ustawienia V/div i time/div (jeśli użytkownik podał)
                ch = self.osc_channel_var.get()
                # zmiana skal V/div
                vs = self.vs_var.get().strip()
                if vs:
                    try:
                        vs_val = float(vs)
                        with self.osc_lock:
                            self.device.oscilloscope.write(f"{ch}:SCALE {vs_val}")
                    except Exception:
                        pass
                # zmiana skali poziomej
                hs = self.hs_var.get().strip()
                unit = self.hs_unit_var.get()
                if hs:
                    try:
                        raw = float(hs)
                        factor = {"s":1.0, "ms":1e-3, "us":1e-6, "ns":1e-9}[unit]
                        val_seconds = raw * factor
                        with self.osc_lock:
                            self.device.oscilloscope.write(f"HOR:MAIN:SCALE {val_seconds}")
                    except Exception:
                        pass

                # pobierz dane z oscyloskopu w locku
                with self.osc_lock:
                    try:
                        self.device.oscilloscope.write(f"DATA:SOU {ch}")
                        self.device.oscilloscope.write("DATA:WIDTH 1")
                        self.device.oscilloscope.write("DATA:ENC RPB")
                        self.device.oscilloscope.write("CURVE?")
                        raw = self.device.oscilloscope.read_raw()
                    except Exception as e:
                        # jeśli I/O padnie, wypisz i zamień raw na None
                        print("Błąd I/O podczas CURVE?:", e)
                        raw = None

                if raw:
                    header_len = 2 + int(raw[1])
                    data = np.frombuffer(raw[header_len:-2], dtype=np.uint8)
                    voltage = (data - self.yoff) * self.ymult + self.yzero
                    time_axis = np.arange(0, len(voltage) * self.xincr, self.xincr)

                    # czasowy
                    self.line_time.set_data(time_axis, voltage)
                    self.ax_time.relim(); self.ax_time.autoscale_view()

                    # FFT
                    N = len(voltage)
                    if N > 0:
                        freq_axis = np.fft.rfftfreq(N, d=self.xincr)
                        fft_mag = np.abs(np.fft.rfft(voltage))
                        self.line_fft.set_data(freq_axis, fft_mag)
                        self.ax_freq.relim(); self.ax_freq.autoscale_view()

                    # draw
                    try:
                        self.canvas.draw()
                    except Exception as e:
                        print("Canvas draw error:", e)
                else:
                    # jeżeli brak raw to rysujemy pusty wykres (nie blokujemy)
                    pass

            except Exception as e:
                print("Błąd przy aktualizacji wykresu:", e)

            time.sleep(PLOT_INTERVAL)

    def on_close(self):
        self._stop_plot_thread = True
        try:
            self.device.close()
        except:
            pass
        self.root.destroy()


# ------------------------- uruchomienie -------------------------
def main():
    root = tk.Tk()
    app = WielkiBratApp(root)
    root.protocol("WM_DELETE_WINDOW", app.on_close)
    root.mainloop()

if __name__ == "__main__":
    main()
