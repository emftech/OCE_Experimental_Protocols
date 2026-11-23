# SBM protocol for the Twin OCT

This folder holds the tools for running the Synchronized BM-mode (SBM) protocol (E. Machado's thesis, Ch. 5) on the Twin OCT. Processing of the acquisitions lives in [`../Codes_SBM_Processing`](../Codes_SBM_Processing). The step-by-step guides mentioned below are kept in the author's lab notes and are not part of this repository.

**Step-by-step guides (Obsidian vault, `Research-GIBIO/LabVIEW/`):**
- **Variant A, external generator (current):** `SBM External Generator Implementation - Step by Step (22-11-2026).md`
- **Variant B, NI DAQ sine:** `Variant_B_DAQ_Sine/SBM_AsyncOCE_DAQ_Excitation_Steps.md` (Spanish) and the 24-10 preparation notes
- **Background and history:** `SBM Protocol Implementation - Twin OCT.md` and `LabView Implementation of SBM Protocol with External Signal Generator Excitation.md`

## Folder layout (22-11-2026)

```
SBM_Protocol/
├── README.md
├── SBM_TwinOCT_Config.m           ← ALL protocol parameters + the VARIANT switch (edit this)
├── SBM_TwinOCT_Waveform.m         ← builds the AO waveform and the timeline (Info)
├── SBM_TwinOCT_LabVIEWNode.m      ← called by the LabVIEW script node; saves the sidecar
├── ScriptNode_GalvoDriver_SBM.txt ← text of the script node (Galvo 3D cluster driver_SBM.vi)
├── Variant_A_ExternalGenerator/   ← Rigol DG922 Pro (1-cycle triggered burst)
│   └── SBM_TwinOCT_GeneratorSettings.m   prints the Rigol + LabVIEW air-puff settings
├── Variant_B_DAQ_Sine/            ← NI DAQ synthesizes the sine on ao3 (AsyncOCE-style)
│   └── SBM_AsyncOCE_DAQ_Excitation_Steps.md
├── Tools/                         ← checks, validation and quick processing
│   ├── SBM_TwinOCT_PlotVectors.m, SBM_TwinOCT_Validate.m, SBM_TwinOCT_Timeline.m
│   ├── SBM_TwinOCT_NullCheck.m, SBM_TwinOCT_CheckExcitation.m, SBM_TwinOCT_ProcessStarter.m
│   ├── SBM_TwinOCT_SmokeTestBin.m, SBM_TwinOCT_SmokeTestNull.m
│   └── SBM_TwinOCT_Recover.m, _Reorganize.m, _Plan.m, _ExportUserWaveform.m, SBM_PhaseCmap.m
└── Validation/                    ← figures written by the validation tools
```

Only the folder itself has to be on the MATLAB path; `SBM_TwinOCT_Config` adds `Tools/` and the variant folders. The LabVIEW script node needs no change: it still does `addpath('<...>\SBM_Protocol')`.

## Choosing the variant (top of the excitation section of `SBM_TwinOCT_Config.m`)

| `Cfg.Variant` | `Cfg.GenTrigger` | What drives the piezo | LabVIEW requirement | AO rows |
|---|---|---|---|---|
| **`"external_generator"`** (default) | **`"puff"`** (default) | Rigol, re-triggered every 57 A-lines by LabVIEW's *Repeated micro puffs* train on ao3 | Enable Puff **ON**, Repeated micro puffs, Puff Frequency 3508.7719 Hz, ON 0.05 ms, Delay 0 (no block-diagram edit) | 3 |
| `"external_generator"` | `"ao3_row"` | Rigol, trigger = TTL row 4 of the SBM waveform | case `"AsyncOCE", "SBM"` in `LabVIEW_Setup_Mirrors_and_Btrig.vi`, Enable Puff OFF | 4 |
| `"external_generator"` | `"counter"` | Rigol, trigger = ÷57 DAQ counter (`SBM_CounterTrigger.vi`) | counter VI running | 3 |
| `"daq_sine"` | — | NI DAQ sine on ao3 → amplifier (no Rigol) | case `"AsyncOCE", "SBM"`, Enable Puff OFF | 4 |

The sidecar stores `Cfg.Variant` and `Cfg.ExcMode`. The processing is identical for all variants.

## Quick start
1. Edit `SBM_TwinOCT_Config.m` (variant, period, A-lines, turn samples, B-scans, frames).
2. Check the vectors with `chk = SBM_TwinOCT_PlotVectors`; it must end with "All pass/fail checks PASSED". Optionally run `SBM_TwinOCT_Validate`.
3. **Variant A:** run `G = SBM_TwinOCT_GeneratorSettings([], [], "DG922Pro")` and set the Rigol and the LabVIEW air-puff panel exactly as printed.
4. In LabVIEW, select the preset `SBM_37A_21B` (Scan Type = SBM, 1197 / 20 / 1 / 57 / 30). The node stops with an error if the GUI and the Config disagree.
5. Before the first acquisition, check on the scope: the trigger period is 285.0 µs, and each trigger gives one sine cycle (Variant A) or a continuous 0–1 V sine on ao3 (Variant B).
6. After each acquisition run `SBM_CheckAcquisitionConsistency(bin)` (in `Codes_SBM_Processing`), then `Main_Processing_SingleExp_SBM_EMF`, which must print EXCITATION DETECTED.
7. Without excitation: set `ExcMode = "none"` (or switch off the excitation) for a null run, then `SBM_TwinOCT_NullCheck(bin)` recommends the B-frame Signal Delay.

## Timing rule
The rule is enforced on the A-line clock, with frame period F = K·P + N_phi. Every frame therefore starts exactly N_phi samples later in the excitation cycle (thesis Eq. 43/47/54). The excitation must be exactly P-periodic on the A-line clock: the generator variants re-trigger a 1-cycle burst every P A-lines, and the DAQ variant is sampled by the same clock. The k-clock cannot be used as a reference.

## Galvo trajectories
Every galvo transition is quintic (C²). Turns longer than 21 samples use stop–dwell–go. The A-line ramps are linear, which keeps the spatial sampling uniform.
