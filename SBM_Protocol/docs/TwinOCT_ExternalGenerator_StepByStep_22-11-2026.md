# SBM External Generator Implementation – Step by Step (22-11-2026)

**Variant A of `SBM_Protocol`:** the Rigol DG922 Pro drives the piezo with a **1-cycle triggered burst**, re-triggered **every 57 A-lines** by LabVIEW's own *Repeated micro puffs* pulse train on **ao3**. This is the same trigger path and cable the QH Multimeridian MB mode already uses (checked on the scope on 20-10), so it needs **no change to the AO channel logic**.

Variant B (NI DAQ synthesizes the sine on ao3) stays in the same folder for later testing; see section "11. Switching to Variant B (DAQ sine) later".


---

## 0. What changed in the code (22-11-2026)

- `SBM_Protocol` has been reorganized. The **variant is chosen at the top of the excitation section of `SBM_TwinOCT_Config.m`**:
  ```matlab
  Cfg.Variant    = "external_generator";   % "external_generator" | "daq_sine"
  Cfg.GenTrigger = "puff";                 % "puff" | "ao3_row" | "counter"
  ```
  The default is `external_generator` + `puff`, which gives `ExcMode = "gen_puff"` and a waveform of **3 AO rows** (X, Y, B-trigger). LabVIEW adds the puff train on ao3.
- Folder layout:
  ```
  SBM_Protocol/
  ├── SBM_TwinOCT_Config.m, SBM_TwinOCT_Waveform.m, SBM_TwinOCT_LabVIEWNode.m, ScriptNode_GalvoDriver_SBM.txt, README.md
  ├── Variant_A_ExternalGenerator/   SBM_TwinOCT_GeneratorSettings.m
  ├── Variant_B_DAQ_Sine/            SBM_AsyncOCE_DAQ_Excitation_Steps.md
  ├── Tools/                         PlotVectors, Validate, Timeline, NullCheck, CheckExcitation, ProcessStarter, SmokeTests, …
  └── Validation/                    figures
  ```
- **The script node text doesn't change.** It still adds only the `SBM_Protocol` folder to the path; the Config adds `Tools/` and the variant folders.
- Verified on 22-11:
  - `SBM_TwinOCT_PlotVectors` passes all checks: 3 rows for Variant A (puff), 4 rows for A (`ao3_row`) and B.
  - The LabVIEW node call with the preset values returns 3 rows, `TotalClusterScans` = 68 229, and saves the sidecar with `ExcMode = gen_puff`.
  - `SBM_TwinOCT_Validate` is unchanged: recovered phase error 6·10⁻¹³ rad, 10.01 m/s for a true 10.00.
  - Both smoke tests pass: synthetic `.bin` at 10.02 m/s; null run recommends delay 33 = the planted lag.

## 1. What you need

- Twin OCT with the master `ATS9360_Master_SSOCT_Trial_3.vi` and the preset `SBM_37A_21B`
- Rigol DG922 Pro, CH1
- Piezo + amplifier, wired as in QH (Rigol CH1 → amplifier → piezo/needle)
- The BNC cable that carries **ao3** to the Rigol's external trigger input (the QH trigger cable)
- Oscilloscope (Rigol MSO5104) with 2 probes, BNC T-pieces
- PC phantom and the thick needle, positioned as for QH
- MATLAB on the Twin PC

## 2. Copy the files to the Twin PC

1. Copy the **whole** `SBM_Protocol` folder (this repository's `SBM_Protocol/`) over the one on the Twin PC.
   - Keep the same location, e.g. `C:\path\to\SBM_Protocol`.
   - Before overwriting, delete the old `.m` files that now live in subfolders. Otherwise old copies at the root shadow the new ones.
2. Also copy `Codes_SBM_Processing` from this repository (it has `SBM_CheckAcquisitionConsistency`).
3. Open `SBM_TwinOCT_Config.m` and check:
   - `Cfg.Variant = "external_generator"`, `Cfg.GenTrigger = "puff"`
   - `Cfg.ExcVpp_V = 1.0`, `Cfg.ExcOffset_V = 0.5` (the QH drive: 1 Vpp, 500 mV)
   - `Cfg.TTL_Width = 10` (puff ON length 0.05 ms), `Cfg.PuffDelay_Alines = 0`
   - protocol: `T_excitation = 285e-6`, `nAscans = 37`, `nTurnSamples = 20`, `nBscans = 21`, `nFrames = 57`, `nPeriods = 4`
4. The script node in `Waveform SubVIs/Galvo 3D cluster driver_SBM.vi` stays as it is. Its two edited paths (the `SBM_Protocol` folder and the sidecar folder) remain valid.

## 3. MATLAB checks on the Twin PC (2 min)

```matlab
cd('C:\path\to\SBM_Protocol')
chk = SBM_TwinOCT_PlotVectors;                          % must end: All pass/fail checks PASSED, Channels: 3
G   = SBM_TwinOCT_GeneratorSettings([], [], "DG922Pro") % prints the Rigol and LabVIEW settings below
```
If `Channels` says 4, the Config isn't on the puff variant (see step 2.3).

## 4. Wiring (identical to QH MB mode)

| From | To |
|---|---|
| DAQ **ao3** BNC (QH trigger) | Rigol **CH1 external trigger input**, the same cable and connector as in QH. **Not** the front *Counter* BNC. |
| Rigol **CH1 output** (G1) | amplifier input → piezo / needle |
| Scope CH1 | ao3, via a T on the trigger cable |
| Scope CH2 | Rigol CH1 output, via a T before the amplifier |

## 5. Program the Rigol DG922 Pro (CH1)

The values match your QH photos of 20-10, except the **period**:

| Step | Setting |
|---|---|
| 1 | **Sine**. Tap the Freq/Period key until it shows **Period** → **282.000 µs** (= 285 µs trigger period − 3 µs dead time to re-arm) |
| 2 | **Ampl 1.000 Vpp**, **Offset 500 mVdc**, Phase **0°**, CH1 Imped **HighZ** |
| 3 | Mode **Burst** → type **NCycle**, **Cycles = 1** (User) |
| 4 | **Source = Ext Leading** (external, rising edge), **Delay 0 s**, burst phase 0° |
| 5 | **Idle Level = 1st Pt** |
| 6 | Burst Period (e.g. 10 ms) is ignored in external-trigger mode; leave it |
| 7 | CH1 **On**. *Store/Recall* → save as `SBM_285us` |

For another excitation period T (frequency sweep), set the Rigol period to **T − 3 µs** (e.g. 220 µs → 217 µs) and the LabVIEW Puff Frequency to **1/T exactly** (e.g. 4545.4545 Hz). `SBM_TwinOCT_GeneratorSettings` prints both for the current Config.

## 6. LabVIEW: air-puff controller and preset

### 6.1 Quick way (no block-diagram edit, enough for today)
1. Run the master and select **Set AcqMode Params = `SBM_37A_21B`**. If it is already selected, pick another preset first and then `SBM_37A_21B`, so the event fires.
2. Check the Scan Parameters:

   | Field | Value |
   |---|---|
   | Scan Type | SBM |
   | A-lines per B-frame | 1197 |
   | B-frames ([3D] Vert) | 1 |
   | No. Cluster | 57 |
   | No. Sync Samples | 20 |
   | B-frame Signal Delay | 30 (or the null-run value, §8) |
   | Scan length H × V | as planned (12 × 12 or 18 × 18 mm) |

3. **After** selecting the preset, set the **Air puff controller** by hand:

   | Control | Value |
   |---|---|
   | **Enable Puff** | **ON** (green) |
   | **Enum** | **Repeated micro puffs** |
   | **Puff Frequency (Hz)** | **3508.7719** (enter every decimal) |
   | **Puff ON Length (ms)** | **0.05** |
   | **Puff Delay (ms)** | **0** |
   | Puff Frequency sweep start / end / step / bin length | 0 |
   | Macro puff voltage | as in QH (it is the TTL level) |
   | Wave Freq / Wave Amp | 0 |

4. Don't touch the preset enum again, or it rewrites Enable Puff OFF.

### 6.2 Permanent way (make the preset write the puff values)
1. Stop the VI. On the front panel, right-click the **Air puff controller** border → **Find → Local Variables**.
2. Step through the results (double-click or Ctrl+G) until the one inside the case **`"SBM_37A_21B"`** of the event `"Set AcqMode Params": Value Change`.
3. Its **Bundle By Name** has *Enable Puff* and *Transducer*. Drag its bottom edge to add **Puff Frequency (Hz)**, **Puff ON Length (ms)** and **Puff Delay (ms)**.
   - Shortcut: copy the Bundle By Name of the `QH` case with its constants, paste it here, and change the values.
4. Constants: Enable Puff **T**; Transducer **Repeated micro puffs**; Puff Frequency **3508.7719** (set the display format to 4 decimals); ON Length **0.05**; Delay **0**.
5. Make sure the Run arrow isn't broken, then **Ctrl+S**.
6. Check the `"SBM"` case that holds a **Diagram Disable** structure (seen on 20-10). If its **Enabled** frame writes the Air puff controller (e.g. Enable Puff = False), delete those writes: they would override the panel on every loop iteration while Scan Type = SBM.

## 7. Scope check before acquiring (mandatory)

Scope: CH1 = ao3, CH2 = Rigol output. Trigger on CH1, rising, 2.5 V, Single. Start a scan with Save off (or a throw-away name).

| Scale | Must see |
|---|---|
| 100 µs/div | CH1 **pulses every 285.0 µs** (*Measure → Period* = 285.0 µs, Freq 3.5088 kHz). CH2 **one 282 µs sine cycle (0–1 V) per pulse**, with ~3 µs flat before the next one |
| 50 ms/div | a **continuous** pulse and sine train for the whole scan (~0.41 s, with the 50-period preamble before the first frame); nothing before or after |

- **Period 280 or 290 µs:** LabVIEW rounded the puff repetition to 56 or 58 A-lines. Set Puff Frequency to **3508.7720** or **3508.7718** and measure again. Keep the value that gives 285.0 µs, and put it in the preset constant.
- **No pulses:** Enable Puff is OFF, or the Enum isn't *Repeated micro puffs*. Check the preset and the Diagram Disable frame (§6.2.6).
- **Pulses but no sine:** check the Rigol (Burst ON, Source Ext, CH1 On) and the trigger cable.
- **Sine cycles missing (gaps):** the Rigol isn't re-arming. Lower its period (e.g. 280 µs) and set `Cfg.Gen_DeadTime_s` to match (5e-6).

## 8. Null run (recommended, once per session)

1. Switch the Rigol CH1 output **OFF**, keeping everything else as is.
2. Acquire one SBM file, e.g. `PC_H_285us_1c_0Vpp_37_20_21_SBM_57_30_12_12_null`.
3. Run `SBM_TwinOCT_NullCheck(bin)`. It reports file and raster OK, frame stability, noise floor, and the **recommended B-frame Signal Delay** (galvo lag from the serpentine shift).
4. If it differs from 30, type the new value into the GUI and use it in the file names (the `Delay` field).
5. Switch the Rigol output back **ON**.

## 9. Acquire

1. File names follow `Experiment_Log_SBM_19_10_2026.xlsx`: `PC_H_285us_1c_1Vpp_37_20_21_SBM_57_30_12_12_r1`, …
2. Start with **one** 285 µs test (r1), and check it (§10) before running the rest.
3. Then tests r2–r5, the B-scan sweep (11 / 41 / 81 B-scans: change `nBscans`, select the matching GUI values that the node asks for) and the frequency sweep (change `T_excitation`, `nTurnSamples`, `nFrames` per the log; Rigol period = T − 3 µs; Puff Frequency = 1/T).
4. After each file, copy the newest `SBM_Waveform_<yyyyMMdd_HHmmss>.mat` from the sidecar folder. The processing matches it to the `.bin` by time stamp, so there's no need to rename it.

## 10. Check and process (Twin PC or Mac)

Folders: `.bin` → `SBM_Experiments/Data`, sidecars → `SBM_Experiments/Sidecars`.

```matlab
[T, ok] = SBM_CheckAcquisitionConsistency(binFile)   % ALL CONSISTENT; ExcMode "gen_puff"; AO rows 3 3
Main_Processing_SingleExp_SBM_EMF                      % set ExpName in section 0
```
The console must show:
- `excitation bin 4: SNR … dB`, clearly above 0
- `detection: … -> EXCITATION DETECTED`

The thick-needle blanking radius is already the default (`Cfg.NeedleRadius_mm = 4.25`). Results: `SBM_Experiments/Results/<name>/` (elastograms, space-time maps, ring profile, cropped and `_Full` videos).

**If NOT DETECTED**:
- Check §7 again, with the scope connected **during** a real acquisition.
- Check that the piezo actually moves the phantom (amplifier on, needle in contact).
- Check the band printout. The real 18/19-10 data showed bulk vibration around 113/178 Hz, which the bulk compensation removes.

## 11. Switching to Variant B (DAQ sine) later

1. In `SBM_TwinOCT_Config.m`, set `Cfg.Variant = "daq_sine"`. This gives 4 AO rows; ao3 carries the 0–1 V sine.
2. In LabVIEW, change the `"AsyncOCE"` case label(s) in `LabVIEW_Setup_Mirrors_and_Btrig.vi` to **`"AsyncOCE", "SBM"`**, and set the preset to Enable Puff **OFF**.
3. Wiring: **ao3 → amplifier input** directly; the Rigol is not used.
4. Full steps: `SBM_Protocol/Variant_B_DAQ_Sine/SBM_AsyncOCE_DAQ_Excitation_Steps.md`.

## 12. Back to QH MB mode

Select the **QH** preset; its case rewrites the puff settings (400 Hz / 0.5 ms / 0.1 ms). Recall the Rigol QH setup (3 cycles). The wiring is the same. `SBM_Protocol` doesn't affect QH.
