# SBM protocol reference

Rules, inputs, outputs, validation tests and lessons learned for the Synchronized BM-mode (SBM) protocol, a novel protocol introduced in E. Machado Ferrer's thesis (Chapter 5) and not yet published in a peer-reviewed journal. For porting to a spectral-domain OCT, see the main [README](../README.md#using-the-code-on-a-spectral-domain-oct).

## The protocol in one page (details: thesis Chapter 5)

SBM acquires the 2D wavefield of a **harmonic** mechanical excitation (needle / piezo, period `T_exc`) with an ordinary **raster BM scan** (N_A A-lines per B-scan, N_B B-scans per frame). It repeats the frame N_F times so that, frame after frame, every pixel is sampled at a slightly later **excitation phase**. Stacking the frames gives, at every pixel, the harmonic sampled in **equivalent time**. One DFT bin per pixel then recovers amplitude and phase. A deterministic per-pixel compensation removes the raster timing (Eq. 59), and 2D phase-gradient estimators give the phase-speed map.

Everything runs on **one clock**: the A-line clock, with period `Ts = 1/line rate`.

| Symbol | Meaning | Rule | Canonical value (200 kHz) |
|---|---|---|---|
| `P = T_exc / Ts` | A-lines per excitation period | **must be an integer** (the excitation is locked to the A-line clock) | 285 µs / 5 µs = **57** |
| `N_A`, `N_turn` | valid A-lines per B-scan, turn samples | `S = N_A + N_turn = k·P` (every B-scan starts at the same phase) | 37 + 20 = 57 |
| `N_B` | B-scans per frame | 11 / 21 / 41 / 81 in the thesis | 21 |
| `N_F`, `N_P` | frames, reconstructed periods | `N_φ = N_P·P/N_F` **integer** and `< P/2` (Eqs. 45–47) | 57, 4 → `N_φ = 4` |
| `F = K·P + N_φ` | frame period (samples) | `F mod P = N_φ`: each frame advances the phase by N_φ samples | 24·57 + 4 = 1372 (6.86 ms) |
| `T_SBM = N_φ·Ts` | equivalent sampling period | the frames sample the harmonic every T_SBM | 20 µs |
| `φ_acq(x,y) = ω(s₁ + τ(x,y))Ts` | acquisition phase of every pixel | known exactly from the timeline; compensated in Eq. 59 | — |

The N_F frames cover exactly N_P periods, so f_exc falls **exactly on DFT bin N_P**: no window and no leakage.

---

## Inputs of the protocol

### `SBM_Protocol/SBM_TwinOCT_Config.m` (edit this file only)

| Field | Unit | Meaning |
|---|---|---|
| `Ts_OCT` | s | **A-line period** = 1/line rate. On SD-OCT: 1 / camera line rate |
| `T_excitation` | s | excitation period. Must be an integer multiple of `Ts_OCT` |
| `nAscans`, `nTurnSamples`, `nBscans` | samples | N_A, N_turn, N_B (raster) |
| `nFrames`, `nPeriods` | — | N_F, N_P (they set N_φ) |
| `horizontalMin_V` … `verticalMax_V` | V | galvo command range (field of view). In the Twin these come from LabVIEW (scan length × V/mm) |
| `MaxAlinesPerBuffer` | A-lines | largest record the acquisition card / frame grabber takes per trigger; sets G, R, C (below) |
| `BscansPerRecord` | — | [] = automatic: G = as many B-scans per record as fit |
| `TrigDelay` | A-lines | **galvo lag**: the B-trigger is sent this many A-lines after the galvo command (measure with `SBM_TwinOCT_NullCheck`) |
| `TrigWidth`, `TrigAmp_V` | A-lines, V | B-trigger pulse |
| `AlazarLatency` | A-lines | delay between the trigger edge and the first recorded A-line (global phase only) |
| `ExcMode` | — | how the excitation is produced (see "Excitation" below) |
| `ExcVpp_V`, `ExcOffset_V` | V | excitation amplitude and offset ("sine": the DAQ output itself) |
| `PreamblePeriods`, `ExcRampPeriods`, `PostPeriods` | periods | excitation on before frame 1 (steady state), fade in/out, tail |
| `FlybackMinSamples`, `TurnReversalSamples`, `RunInSamples` | samples | galvo smoothness (quintic C² transitions; stop-dwell-go for long turns) |
| `TTL_*`, `Gen_*`, `PuffDelay_Alines` | — | only for the external-generator modes |

The generator **refuses** infeasible settings with an explicit error: non-integer P or N_φ, N_φ ≥ P/2, a record larger than `MaxAlinesPerBuffer`, or a waveform beyond ±10 V.

### Excitation

The excitation variant is selected at the top of the excitation section of `SBM_TwinOCT_Config.m` with `Cfg.Variant` (`"external_generator"`, default, or `"daq_sine"`) and `Cfg.GenTrigger` (`"puff"`, default, `"ao3_row"` or `"counter"`); `ExcMode` is derived from them. See the variant table in [`SBM_Protocol/README.md`](../SBM_Protocol/README.md). `ExcMode = "none"` gives a null run (no excitation).

### Feasible settings at other line rates (computed with this generator; N_A = 37, N_F = P, N_P = 4, turn ≥ 100 µs)

| Line rate | Target f | P | f_exc | N_turn | T_SBM | Acq. window (21 B-scans) |
|---|---|---|---|---|---|---|
| 20 kHz | 2 kHz | 10 | 2000 Hz | 3 | 200 µs | 437 ms |
| 50 kHz | 3.5 kHz | 14 | 3571 Hz | 5 | 80 µs | 260 ms |
| 76 kHz | 3.5 kHz | 22 | 3455 Hz | 29 | 53 µs | 421 ms |
| 100 kHz | 3.5 kHz | 29 | 3448 Hz | 21 | 40 µs | 379 ms |
| 147 kHz | 3.5 kHz | 42 | 3500 Hz | 47 | 27 µs | 529 ms |
| 200 kHz | 3.5 kHz | 57 | 3509 Hz | 20 | 20 µs | 390 ms |
| 250 kHz | 3.5 kHz | 71 | 3521 Hz | 34 | 16 µs | 464 ms |

At low line rates, keep P ≥ 9 so that N_φ = 4 < P/2. N_A can be changed, but S = N_A + N_turn must remain a multiple of P.

---

## Outputs of the protocol generator

`[waveform, Info] = SBM_TwinOCT_Waveform(Cfg)`

- **`waveform`** [rows × Ntot], volts, one column per A-line clock tick. Output it with the AO **sample clock = A-line clock**.
  - row 1: galvo X
  - row 2: galvo Y
  - row 3: **B-trigger** (5 V pulses that start each record)
  - row 4: excitation (only for "sine" / "gen_ao3" / "ttl")
- **`Info`**: the exact timeline, saved as the **sidecar** `.mat` next to every acquisition. The processing uses it.

| Field | Meaning |
|---|---|
| `P, Nphi, S, F, K, G, R, C, Nbody, Ntot` | protocol integers (R = A-lines per record, C = records per frame) |
| `FrameStart`, `TrigIdx` [C × N_F] | sample index of each frame start and each record trigger |
| `PixelIdx` [N_B × N_A] | raster pixel → A-line index inside the frame record (serpentine handled) |
| `TauSamples`, **`PhiAcq`** [N_B × N_A] | acquisition delay and phase of every pixel (Eq. 57) |
| `FramePhaseSamples`, `Ts_SBM`, `f_exc`, `AcqWindow_s`, `Duration_s` | equivalent-time bookkeeping |
| `xAxis_V`, `yAxis_V`, `Overshoot_V`, `RunInSamples`, `FlybackSamples` | geometry and galvo information |
| `LabVIEW` | the values that the acquisition software must be set to (record length, records per frame, number of frames, sync samples, trigger delay) |

`T = SBM_TwinOCT_Timeline(Period_us, NA, NB, NP, NF, FOV_mm)` returns the same information per recorded A-line (`T.n`, `T.t_s`, `T.x_mm`, `T.y_mm`) and per pixel (`T.Pix.t_s` [N_B × N_A × N_F]). Use it to drive simulations with the real hardware timing.

**Acquired data expected by the processing:** N_F frames, each of R·C consecutive A-lines (C records of R A-lines started by the B-triggers), stored in acquisition order, plus the sidecar.

## Validation and tests (run these first on the new system)

| Test | Command | Expected (11-11-2026 run) |
|---|---|---|
| Waveform checks | `chk = SBM_TwinOCT_PlotVectors` | "All pass/fail checks PASSED" |
| Recovery from the timeline | `run SBM_TwinOCT_Validate.m` | phase error 1.5 rad raw → ~1e-13 rad recovered; speed 10.01 m/s (true 10.00) |
| Synthetic raw file | `run SBM_TwinOCT_SmokeTestBin.m` | phase error 0.009 rad, 10.02 m/s |
| Null run / galvo lag | `run SBM_TwinOCT_SmokeTestNull.m` | recommended delay 33 = planted lag |
| Full processing chain | `Codes_SBM_Processing/Tests/SBM_Test_Pipeline.m` (7 synthetic acquisitions, ~5 min) | 285 µs, 21 B-scans: 2D PG-TV + DF ≈ 10.61 m/s (true 10.54), Ep ≈ 1.5 %, every run "EXCITATION DETECTED", needle at (0.70, 1.10) mm |

**Experimental sequence that worked best on the Twin OCT:**
1. `SBM_TwinOCT_PlotVectors` → PASS.
2. **Scope** on the excitation channel during a scan: a sine of period P·Ts, present for the whole acquisition.
3. **Null run** (`ExcMode "none"`) → `SBM_TwinOCT_NullCheck` gives the galvo lag (set `TrigDelay`), the noise floor and the frame stability.
4. First excitation run → `SBM_CheckAcquisitionConsistency` (header ↔ sidecar ↔ name), then the single-experiment script. The console must print **`EXCITATION DETECTED`**.

---

## Lessons learned on the Twin OCT (Oct–Nov 2026)

- **First phantom tests had no excitation.** The acquisition software used a 3-channel AO task and silently dropped the 4th (excitation) row. The detection test in `SBM_HarmonicRecovery` flagged it (bin-N_P level at the noise floor, no spatial coherence). Always check the excitation channel on a scope during a real scan.
- **Bulk motion.** Frames are milliseconds apart (6.86 ms at 200 kHz), so slow table or phantom vibration gives lag-1 phase changes of ±π that alias over the whole frame spectrum. `SBM_BulkCompensation` (circular moving mean over one raster line) brings this from 1.8 rad to 0.2 rad. It is enabled by default.
- **Sensitivity.** With 57 frames and about 0.2 rad of noise per frame pair, the detection floor was about 30–40 nm of displacement at the source. Use enough drive amplitude, more frames (larger N_P), or coherent averaging of repeated runs (the timeline is deterministic).
- **Thick needle near field.** Simulations show 2D PG is biased within r ≈ 4.25 mm of a 3.5 mm-wide needle (`Cfg.NeedleRadius_mm`). Metrics use r ≥ R + 1 mm (far-field ROI).
- **Only harmonics of f_exc are physical in equivalent time.** Any non-locked vibration aliases through the frame rate onto arbitrary bins, so a peak away from bin N_P is not the wave.

---

