"""
generator/generacja.py

Moduł odpowiadający za:
 - połączenie z generatorem sygnałów (Siglent/zgodny SCPI),
 - podstawowe funkcje (SINE, SQUARE, RAMP, PULSE, NOISE, ARB, DC),
 - ustawienia wyjścia (LOAD, WAVE_COMBINE),
 - funkcja manual_sweep która wykonuje ręczny sweep częstotliwości,
   odczytuje dane z podanego oscyloskopu (przekazanego jako obiekt pyvisa)
   i zapisuje statystyki (Vpp, Mean, Min, Max, StdDev) do CSV.
"""

import time
import csv
import numpy as np
import pyvisa
import threading

# opcjonalny lock (można podać z main.py lub użyć lokalnego)
# jeśli chcesz, możesz nadpisać osc_lock z main.py (np. pass w manual_sweep)
osc_lock = threading.Lock()


class GeneratorDevice:
    def __init__(self, generator_address=None, oscilloscope_address=None, timeout=10000):
        """
        Nie otwieramy połączeń automatycznie, dopiero przez connect_*(),
        żeby import modułu nie blokował programu jeśli brak urządzeń.
        """
        self.rm = pyvisa.ResourceManager()
        self.generator_address = generator_address
        self.oscilloscope_address = oscilloscope_address
        self.generator = None
        self.oscilloscope = None
        self.timeout = timeout

    # ---------- połączenia ----------
    def connect_generator(self, address=None):
        addr = address or self.generator_address
        if not addr:
            raise ValueError("Brak adresu generatora.")
        self.generator = self.rm.open_resource(addr)
        self.generator.timeout = self.timeout
        return self.generator.query("*IDN?").strip()

    def connect_oscilloscope(self, address=None):
        addr = address or self.oscilloscope_address
        if not addr:
            raise ValueError("Brak adresu oscyloskopu.")
        self.oscilloscope = self.rm.open_resource(addr)
        self.oscilloscope.timeout = self.timeout
        return self.oscilloscope.query("*IDN?").strip()

    # ---------- generator: podstawowe fale ----------
    def SINE(self, channel="C1", frequency=1.0, amplitude=1.0, offset=0.0, phase=0.0):
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:BSWV WVTP,SINE")
        self.generator.write(f"{channel}:BSWV FRQ,{frequency}")
        self.generator.write(f"{channel}:BSWV AMP,{amplitude}")
        self.generator.write(f"{channel}:BSWV OFST,{offset}")
        self.generator.write(f"{channel}:BSWV PHSE,{phase}")
        self.generator.write(f"{channel}:OUTP ON")

    def SQUARE(self, channel="C1", frequency=1.0, amplitude=1.0, offset=0.0, phase=0.0, duty=50.0):
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:BSWV WVTP,SQUARE")
        self.generator.write(f"{channel}:BSWV FRQ,{frequency}")
        self.generator.write(f"{channel}:BSWV AMP,{amplitude}")
        self.generator.write(f"{channel}:BSWV OFST,{offset}")
        self.generator.write(f"{channel}:BSWV PHSE,{phase}")
        self.generator.write(f"{channel}:BSWV DUTY,{duty}")
        self.generator.write(f"{channel}:OUTP ON")

    def RAMP(self, channel="C1", frequency=1.0, amplitude=1.0, offset=0.0, phase=0.0, symmetry=50.0):
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:BSWV WVTP,RAMP")
        self.generator.write(f"{channel}:BSWV FRQ,{frequency}")
        self.generator.write(f"{channel}:BSWV AMP,{amplitude}")
        self.generator.write(f"{channel}:BSWV OFST,{offset}")
        self.generator.write(f"{channel}:BSWV PHSE,{phase}")
        self.generator.write(f"{channel}:BSWV SYM,{symmetry}")
        self.generator.write(f"{channel}:OUTP ON")

    def PULSE(self, channel="C1", frequency=1.0, duty=50.0, amplitude=1.0, offset=0.0,
              pulse_width=0.001, rise_edge=1e-6, fall=1e-6, delay=0.0):
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:BSWV WVTP,PULSE")
        self.generator.write(f"{channel}:BSWV FRQ,{frequency}")
        self.generator.write(f"{channel}:BSWV DUTY,{duty}")
        self.generator.write(f"{channel}:BSWV AMP,{amplitude}")
        self.generator.write(f"{channel}:BSWV OFST,{offset}")
        self.generator.write(f"{channel}:BSWV WIDTH,{pulse_width}")
        self.generator.write(f"{channel}:BSWV RISE,{rise_edge}")
        self.generator.write(f"{channel}:BSWV FALL,{fall}")
        self.generator.write(f"{channel}:BSWV DLY,{delay}")
        self.generator.write(f"{channel}:OUTP ON")

    def NOISE(self, channel="C1", stdev=0.5, mean=0.0):
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:BSWV WVTP,NOISE")
        self.generator.write(f"{channel}:BSWV STDEV,{stdev}")
        self.generator.write(f"{channel}:BSWV MEAN,{mean}")
        self.generator.write(f"{channel}:OUTP ON")

    def ARB(self, channel="C1", frequency=1.0, amplitude=1.0, offset=0.0, phase=0.0):
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:BSWV WVTP,ARB")
        self.generator.write(f"{channel}:BSWV FRQ,{frequency}")
        self.generator.write(f"{channel}:BSWV AMP,{amplitude}")
        self.generator.write(f"{channel}:BSWV OFST,{offset}")
        self.generator.write(f"{channel}:BSWV PHSE,{phase}")
        self.generator.write(f"{channel}:OUTP ON")

    def DC(self, channel="C1", offset=0.0):
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:BSWV WVTP,DC")
        self.generator.write(f"{channel}:BSWV OFST,{offset}")
        self.generator.write(f"{channel}:OUTP ON")

    # ---------- wyjście ----------
    def WAVE_COMBINE(self, channel="C1", state=True):
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:CMBN STATE,{'ON' if state else 'OFF'}")

    def SET_LOAD(self, channel="C1", load='50'):
        """
        LOAD: '50', 'HZ' (High-Z), albo numeric 50..100000
        """
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        if isinstance(load, (int, float)):
            if 50 <= load <= 100000:
                self.generator.write(f"{channel}:OUTP LOAD,{load}")
            else:
                raise ValueError("LOAD musi być 50..100000")
        elif isinstance(load, str) and load.upper() in ("50", "HZ"):
            self.generator.write(f"{channel}:OUTP LOAD,{load.upper()}")
        else:
            raise ValueError("LOAD musi być '50', 'HZ' lub liczba 50..100000")

    def SWEEP_hw(self, channel, start_freq, stop_freq, sweep_time, mode, direction, wave_type, amplitude, offset, frequency):
        """
        Ustawienia sprzętowego sweep w generatorze (opcjonalne, nie używane w manual_sweep).
        """
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        self.generator.write(f"{channel}:BSWV WVTP,{wave_type}")
        self.generator.write(f"{channel}:BSWV FRQ,{frequency}")
        self.generator.write(f"{channel}:BSWV AMP,{amplitude}")
        self.generator.write(f"{channel}:BSWV OFST,{offset}")
        self.generator.write(f"{channel}:SWWV STATE,ON")
        self.generator.write(f"{channel}:SWWV TIME,{sweep_time}")
        self.generator.write(f"{channel}:SWWV START,{start_freq}")
        self.generator.write(f"{channel}:SWWV STOP,{stop_freq}")
        self.generator.write(f"{channel}:SWWV SWMD,{mode}")
        self.generator.write(f"{channel}:SWWV DIR,{direction}")
        self.generator.write(f"{channel}:SWWV TRSR,INT")
        self.generator.write(f"{channel}:OUTP ON")

    # ---------- manual_sweep ----------
    def manual_sweep(self,
                     channel,
                     start_freq,
                     stop_freq,
                     step,
                     amplitude,
                     offset,
                     delay,
                     oscilloscope=None,
                     ymult=None,
                     yzero=None,
                     yoff=None,
                     xincr=None,
                     osc_channel="CH1",
                     csv_file="manual_sweep_output.csv",
                     wave_type="SINE",
                     osc_lock_local=None,
                     progress_callback=None):
        """
        Ręczny sweep: iterujemy po częstotliwości od start_freq do stop_freq
        co 'step' i dla każdej częstotliwości:
          - ustawiamy generator,
          - czekamy `delay` sec,
          - pobieramy przebieg z oscyloskopu (raw CURVE) i liczymy statystyki:
            Vpp, Mean, Min, Max, StdDev,
          - zapisujemy do CSV.
        """

        # walidacja
        if self.generator is None:
            raise RuntimeError("Generator niepodłączony.")
        if oscilloscope is None:
            raise ValueError("Trzeba podać obiekt oscilloscope (pyvisa.Resource).")
        # jeśli brak parametrów skalowania (ymult itp.), spróbujemy je pobrać z oscilloscope
        if ymult is None:
            ymult = float(oscilloscope.query("WFMPRE:YMULT?"))
        if yzero is None:
            yzero = float(oscilloscope.query("WFMPRE:YZERO?"))
        if yoff is None:
            yoff = float(oscilloscope.query("WFMPRE:YOFF?"))
        if xincr is None:
            xincr = float(oscilloscope.query("WFMPRE:XINCR?"))

        # wybierz lock (użyj przekazanego lub globalnego)
        lock = osc_lock_local if osc_lock_local is not None else osc_lock

        with open(csv_file, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Frequency_Hz", "Vpp_V", "Mean_V", "Min_V", "Max_V", "StdDev_V"])

            # upewniamy się, że to floaty
            start = float(start_freq)
            stop = float(stop_freq)
            step = float(step)
            freq = start
            step_count = 0

            while freq <= stop + 1e-12:
                # ustaw generator
                self.generator.write(f"{channel}:BSWV WVTP,{wave_type}")
                self.generator.write(f"{channel}:BSWV FRQ,{freq}")
                self.generator.write(f"{channel}:BSWV AMP,{amplitude}")
                self.generator.write(f"{channel}:BSWV OFST,{offset}")
                self.generator.write(f"{channel}:OUTP ON")

                # daj czas na ustabilizowanie
                time.sleep(delay)

                # odczyt surowy z oscyloskopu w bloku lock
                with lock:
                    try:
                        oscilloscope.write(f"DATA:SOU {osc_channel}")
                        oscilloscope.write("DATA:WIDTH 1")
                        oscilloscope.write("DATA:ENC RPB")
                        oscilloscope.write("CURVE?")
                        raw = oscilloscope.read_raw()
                    except Exception as e:
                        # jeśli CURVE? nie działa, spróbuj pomiarów sprzętowych (MEAS?)
                        raw = None
                        print("Błąd przy odczycie CURVE?:", e)

                # parsowanie - jeśli raw OK to przetwarzamy; jeśli nie, próbujemy komend MEAS
                if raw:
                    header_len = 2 + int(raw[1])
                    data = np.frombuffer(raw[header_len:-2], dtype=np.uint8)
                    voltage = (data - yoff) * ymult + yzero
                    v_max = float(np.max(voltage))
                    v_min = float(np.min(voltage))
                    v_pp = v_max - v_min
                    v_mean = float(np.mean(voltage))
                    v_std = float(np.std(voltage))
                else:
                    # fallback na MEAS: jeśli dostępne na oscyloskopie
                    try:
                        # niektóre oscyloskopy potrafią zwrócić wartości MEAS:VPP? etc.
                        with lock:
                            v_pp = float(oscilloscope.query("MEAS:VPP?"))
                            v_mean = float(oscilloscope.query("MEAS:MEAN?"))
                            v_min = float(oscilloscope.query("MEAS:MIN?"))
                            v_max = float(oscilloscope.query("MEAS:MAX?"))
                            v_std = float(oscilloscope.query("MEAS:STDDEV?"))
                    except Exception as e:
                        print("Fallback MEAS też nie zadziałał:", e)
                        v_pp = v_mean = v_min = v_max = v_std = float('nan')

                writer.writerow([freq, v_pp, v_mean, v_min, v_max, v_std])
                step_count += 1
                if progress_callback:
                    progress_callback(step_count, freq, v_pp, v_mean, v_min, v_max, v_std)

                freq += step

    # ---------- zamknięcie ----------
    def close(self):
        try:
            if self.generator:
                self.generator.close()
        except:
            pass
        try:
            if self.oscilloscope:
                self.oscilloscope.close()
        except:
            pass


# helper init (łatwiej importować z main.py)
def init_device(generator_address=None, oscilloscope_address=None, timeout=10000):
    dev = GeneratorDevice(generator_address, oscilloscope_address, timeout=timeout)
    return dev
