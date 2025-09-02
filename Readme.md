# Bulbulator – Generator & Oscilloscope GUI

A desktop application written in **Python (Tkinter + Matplotlib + PyVISA)** for controlling a **Siglent SDG signal generator** and a **oscilloscope**.  
Project developed as part of internship / student practice.

---

## ✨ Features

- **Full GUI (Tkinter)**
  - Left panel: generator control (waveforms, parameters, sweep, manual sweep, protection, impedance, IP settings).
  - Right panel: real-time plots (live waveform + FFT).
- **Waveform control**
  - SINE, SQUARE, RAMP, PULSE, NOISE, ARB, DC.
  - Parameters: frequency, amplitude, offset, duty cycle, etc.
- **Sweep modes**
  - Automatic sweep (`:SWWV`).
  - Manual sweep with CSV logging (Vpp, Mean, Min, Max, StdDev).
- **Oscilloscope integration**
  - Fetch waveform via SCPI (`CURVE?`).
  - Live plot in Matplotlib (embedded in Tkinter).
  - Automatic axis scaling.
  - Measurement queries (`MEAS:VPP?`, `MEAS:MEAN?`, `MEAS:MIN?`, `MEAS:MAX?`, `MEAS:STDDEV?`).
- **Extras**
  - OVP (Over-Voltage Protection).
  - Impedance configuration (`OUTP LOAD`).
  - LAN/IP remote settings (`SYST:COMM:LAN:IPADDR`).
  - Dual-channel / Wave Combine mode.

---

## 📦 Requirements

- Python **3.10+**  
- Dependencies (install via `pip install -r requirements.txt`):
