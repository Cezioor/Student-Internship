import pyvisa
import time

# --- KONFIG: adres generatora (Siglent SDG) ---
GENERATOR_ADDRESS = "TCPIP::192.168.0.200::INSTR"

rm = pyvisa.ResourceManager()
device = rm.open_resource(GENERATOR_ADDRESS)
device.timeout = 10000

try:
    print("GEN *IDN?:", device.query("*IDN?").strip())
except Exception as e:
    print("Nie udało się odczytać *IDN? z generatora:", e)


# ========= podstawowe fale =========

def SINE(channel: str = 'C1', frequency: float = 1.0, amplitude: float = 1.0, offset: float = 0.0, phase: float = 0.0):
    device.write(f"{channel}:BSWV WVTP,SINE")
    device.write(f"{channel}:BSWV FRQ,{frequency}")
    device.write(f"{channel}:BSWV AMP,{amplitude}")
    device.write(f"{channel}:BSWV OFST,{offset}")
    device.write(f"{channel}:BSWV PHSE,{phase}")
    device.write(f"{channel}:OUTP ON")

def SQUARE(channel: str = 'C1', frequency: float = 1.0, amplitude: float = 1.0, offset: float = 0.0, phase: float = 0.0, duty: float = 50.0):
    device.write(f"{channel}:BSWV WVTP,SQUARE")
    device.write(f"{channel}:BSWV FRQ,{frequency}")
    device.write(f"{channel}:BSWV AMP,{amplitude}")
    device.write(f"{channel}:BSWV OFST,{offset}")
    device.write(f"{channel}:BSWV PHSE,{phase}")
    device.write(f"{channel}:BSWV DUTY,{duty}")
    device.write(f"{channel}:OUTP ON")

def RAMP(channel: str = 'C1', frequency: float = 1.0, amplitude: float = 1.0, offset: float = 0.0, phase: float = 0.0, symmetry: float = 50.0):
    device.write(f"{channel}:BSWV WVTP,RAMP")
    device.write(f"{channel}:BSWV FRQ,{frequency}")
    device.write(f"{channel}:BSWV AMP,{amplitude}")
    device.write(f"{channel}:BSWV OFST,{offset}")
    device.write(f"{channel}:BSWV PHSE,{phase}")
    device.write(f"{channel}:BSWV SYM,{symmetry}")
    device.write(f"{channel}:OUTP ON")

def PULSE(channel: str = 'C1', frequency: float = 1.0, duty: float = 50.0, amplitude: float = 1.0, offset: float = 0.0,
          pulse_width: float = 0.001, rise_edge: float = 0.000001, fall: float = 0.000001, delay: float = 0.0):
    device.write(f"{channel}:BSWV WVTP,PULSE")
    device.write(f"{channel}:BSWV FRQ,{frequency}")
    device.write(f"{channel}:BSWV DUTY,{duty}")
    device.write(f"{channel}:BSWV AMP,{amplitude}")
    device.write(f"{channel}:BSWV OFST,{offset}")
    device.write(f"{channel}:BSWV WIDTH,{pulse_width}")
    device.write(f"{channel}:BSWV RISE,{rise_edge}")
    device.write(f"{channel}:BSWV FALL,{fall}")
    device.write(f"{channel}:BSWV DLY,{delay}")
    device.write(f"{channel}:OUTP ON")

def NOISE(channel: str = 'C1', stdev: float = 0.5, mean: float = 0.0):
    device.write(f"{channel}:BSWV WVTP,NOISE")
    device.write(f"{channel}:BSWV STDEV,{stdev}")
    device.write(f"{channel}:BSWV MEAN,{mean}")
    device.write(f"{channel}:OUTP ON")

def ARB(channel: str = 'C1', frequency: float = 1.0, amplitude: float = 1.0, offset: float = 0.0, phase: float = 0.0):
    device.write(f"{channel}:BSWV WVTP,ARB")
    device.write(f"{channel}:BSWV FRQ,{frequency}")
    device.write(f"{channel}:BSWV AMP,{amplitude}")
    device.write(f"{channel}:BSWV OFST,{offset}")
    device.write(f"{channel}:BSWV PHSE,{phase}")
    device.write(f"{channel}:OUTP ON")

def DC(channel: str = 'C1', offset: float = 0.0):
    device.write(f"{channel}:BSWV WVTP,DC")
    device.write(f"{channel}:BSWV OFST,{offset}")
    device.write(f"{channel}:OUTP ON")


# ========= opcje wyjścia =========

def WAVE_COMBINE(channel='C1', state=True):
    device.write(f"{channel}:CMBN STATE,{'ON' if state else 'OFF'}")

def SET_LOAD(channel='C1', load='50'):
    """
    Ustawienie impedancji wyjścia generatora: '50', 'HZ' (High Z) lub wartość 50–100000 Ω.
    """
    if isinstance(load, (int, float)):
        if 50 <= load <= 100000:
            device.write(f"{channel}:OUTP LOAD,{load}")
        else:
            raise ValueError("Impedancja poza zakresem 50–100000 Ω.")
    elif isinstance(load, str) and load.upper() in ('50', 'HZ'):
        device.write(f"{channel}:OUTP LOAD,{load.upper()}")
    else:
        raise ValueError("LOAD musi być '50', 'HZ' albo liczbą 50–100000.")


# ========= sweep (ustawienia sprzętowego sweepa) =========

def SWEEP(channel, start_freq, stop_freq, sweep_time, mode, direction, wave_type, amplitude, offset, frequency):
    """
    Ustawia tryb SWWV w generatorze. (Używane sporadycznie – manual_sweep realizujemy pętlą po FRQ.)
    """
    device.write(f"{channel}:BSWV WVTP,{wave_type}")
    device.write(f"{channel}:BSWV FRQ,{frequency}")
    device.write(f"{channel}:BSWV AMP,{amplitude}")
    device.write(f"{channel}:BSWV OFST,{offset}")
    device.write(f"{channel}:SWWV STATE,ON")
    device.write(f"{channel}:SWWV TIME,{sweep_time}")
    device.write(f"{channel}:SWWV START,{start_freq}")
    device.write(f"{channel}:SWWV STOP,{stop_freq}")
    device.write(f"{channel}:SWWV SWMD,{mode}")
    device.write(f"{channel}:SWWV DIR,{direction}")
    device.write(f"{channel}:SWWV TRSR,INT")
    device.write(f"{channel}:OUTP ON")


# ========= manual_sweep (pętla po FRQ + pomiar na oscyloskopie) =========

def manual_sweep(
    channel: str,
    start_freq: float,
    stop_freq: float,
    step: float,
    amplitude: float,
    offset: float,
    delay: float,
    oscilloscope,          # uchwyt VISA do oscyloskopu (z main.py)
    ymult: float,
    yzero: float,
    yoff: float,
    xincr: float,
    osc_channel: str = "CH1",
    csv_file: str = "manual_sweep_output.csv",
    wave_type: str = "SINE",
    osc_lock=None,         # opcjonalnie threading.Lock() z main.py
    progress_callback=None
):
    """
    Ustawia generator na kolejne częstotliwości i dla każdej pobiera przebieg z oscyloskopu,
    liczy statystyki (Vpp, mean, min, max, std) i zapisuje do CSV.
    """
    import csv
    import numpy as np

    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Frequency_Hz", "Vpp_V", "Mean_V", "Min_V", "Max_V", "StdDev_V"])

        freq = float(start_freq)
        while freq <= float(stop_freq) + 1e-12:
            # ustaw generator (bez locka – to inny instrument)
            device.write(f"{channel}:BSWV WVTP,{wave_type}")
            device.write(f"{channel}:BSWV FRQ,{freq}")
            device.write(f"{channel}:BSWV AMP,{amplitude}")
            device.write(f"{channel}:BSWV OFST,{offset}")
            device.write(f"{channel}:OUTP ON")

            time.sleep(delay)

            # pobierz dane z oscyloskopu (w locku jeśli przekazano)
            if osc_lock is None:
                _ctx = _DummyCtx()
            else:
                _ctx = osc_lock

            with _ctx:
                oscilloscope.write(f"DATA:SOU {osc_channel}")
                oscilloscope.write("DATA:WIDTH 1")
                oscilloscope.write("DATA:ENC RPB")
                oscilloscope.write("CURVE?")
                raw = oscilloscope.read_raw()

            # dekoduj
            header_len = 2 + int(raw[1])
            data = np.frombuffer(raw[header_len:-2], dtype=np.uint8)
            voltage = (data - yoff) * ymult + yzero

            v_max  = float(np.max(voltage))
            v_min  = float(np.min(voltage))
            v_pp   = v_max - v_min
            v_mean = float(np.mean(voltage))
            v_std  = float(np.std(voltage))

            writer.writerow([freq, v_pp, v_mean, v_min, v_max, v_std])

            if progress_callback:
                progress_callback(freq, v_pp)

            freq += float(step)


class _DummyCtx:
    """Prosty kontekst, gdy nie przekazano locka."""
    def __enter__(self): return self
    def __exit__(self, exc_type, exc, tb): return False
