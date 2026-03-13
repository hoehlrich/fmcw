# FMCW Radar Project Overview

## Goal
Simulate and implement an FMCW radar signal processing pipeline, culminating in
real hardware validation.

## Background
FMCW radar works by transmitting a frequency-swept chirp, mixing the delayed
return with the live TX signal, and analyzing the resulting IF signal. The IF
frequency encodes target range; Doppler shift across multiple chirps encodes
velocity.

---

## Phase 1 — Range Simulation (MATLAB)
**Goal:** Simulate a single chirp, a target return, and recover target range via FFT.

1. Define radar parameters: bandwidth B, chirp duration T, carrier frequency f0, sample rate fs
2. Generate time vector `t` and TX chirp signal
3. Simulate target return as time-delayed chirp (delay tao = 2R/c)
4. Mix TX and RX signals (elementwise multiply) to get IF signal
5. FFT the IF signal and identify range peak
6. Verify: IF frequency should equal (B/T) × τ
7. Add a second target and verify range resolution

## Phase 2 — Range-Doppler Map (MATLAB)
**Goal:** Simulate moving targets and recover velocity.

1. Simulate a sequence of chirps with a moving target
2. Apply range FFT to each chirp
3. Apply Doppler FFT across chirps (2D FFT)
4. Generate and interpret range-Doppler map
5. Verify recovered velocity matches simulated target motion

## Phase 3 — Hardware Validation
**Goal:** Run the processing pipeline on real radar data.

1. Borrow FMCW hardware (via ISS professor)
2. Capture IF signal from known stationary target
3. Run Phase 1 pipeline on real data, validate range
4. Extend to moving target and validate Doppler

## Notes
- Reach out to professor *after* Phase 1 is working — gives you something concrete to show
- MATLAB surface area needed: `fft()`, `fftshift()`, `plot()`, basic complex arithmetic, elementwise ops
- Prior work (IIR filter, spectrogram) maps directly — same discrete vector / frequency domain mental model
