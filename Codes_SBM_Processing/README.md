# Codes_SBM_Processing

Processing of the **Synchronized BM-mode (SBM)** experiments of the Twin OCT. It goes from the raw LabVIEW `.bin` (+ its sidecar) to the 2D PG / 2D PG-TV elastograms, without and with directional filtering. It is built like the Multimeridian MB-mode pipeline (`Main_Processing_*_MultimerMB_EMF.m`, `Codes/MultimerMB_Pipeline`) and follows the thesis SBM chain (Chapter 5, `Sampling_Acq_Sim_SBM_v4`).

Author: Estéfano Machado Ferrer (with Claude Code), 19-10-2026.

## Folder layout expected

```
OCE_Experimental_Protocols/   (repository root)
├── Codes/                       (MB codes, reused: FindSurface, Loupas2D_Fast, OCE_PG2D, TV, MMB_* plots/videos)
├── Codes_SBM_Processing/        (this folder)
└── SBM_Experiments/
    ├── Data/                    *.bin (names from Experiment_Log_SBM_*.xlsx) [+ <name>_Sidecar.mat]
    ├── Sidecars/                SBM_Waveform_yyyyMMdd_HHmmss.mat / SBM_Waveform_LAST.mat (LabVIEW node)
    └── Results/                 <name>/ ... and Summary_SBM/
```

## Run

- **One experiment:** set `ExpName` in section 0 of `Main_Processing_SingleExp_SBM_EMF.m` and run it.
- **Everything:** run `Main_Processing_Batch_SBM_EMF`. Processed experiments are skipped unless `Overwrite = true`. Use `SummaryOnly = true` to rebuild the summary only, and `WorkerID`/`NumWorkers` to split the list between several MATLAB processes.
- **Test:** `Tests/SBM_Test_Pipeline` builds 7 synthetic `.bin` files (known c = 10.54 m/s, needle at (0.70, 1.10) mm, 30 % boundary reflection) and runs the batch on them.

## MB step → SBM step

| MB (`Main_Processing_SingleExp_MultimerMB_EMF`) | SBM (`Main_Processing_SingleExp_SBM_EMF`) |
|---|---|
| header, memmap, DC spectrum, B-mode, `FindSurface` | same (frame-averaged B-mode along the raster path) + a guard against noise picks |
| Loupas lag-1 along the A-lines of an M-scan | Loupas lag-1 **across frames** (equivalent time T_SBM = N_φ·Ts), circular |
| narrowband Gaussian filter around f_exc | **exact DFT bin N_P** (the N_F frames span exactly N_P periods) |
| (none) | **acquisition-phase compensation, Eq. 59**, with the exact φ_acq from the sidecar |
| needle fit along the meridians | coherence fit on the 2D recovered phasor (`SBM_EstimateSource`) |
| outward k-f filter per meridian | 2D sector-wise outward k-space filter (`SBM_DirectionalFilterOutward2D`) |
| interpolation of the sparse meridians | spline onto a square ≥ 100 × 100 grid (`SBM_InterpolateSquare` = `Prepare_SBM_2DPhaseGradient`) |
| `OCE_PG_ProcessMB` (1D/2D PG, PG-TV) | `OCE_PG2D` (2D PG, 2D PG-TV), as in the thesis SBM results |
| grids, space-time maps, depth/en face videos, `Results_MMB.mat` | same, `Results_SBM.mat` |
| batch: best rep, 4×4 grids, metrics vs meridians, dispersion | batch: test repeatability, B-scan sweep grids, metrics vs N_B vs MB, dispersion vs MB, session needle |

## Functions (`SBM_Pipeline/`)

| Function | Purpose |
|---|---|
| `SBM_ParseExperimentName` | `PC_H_285us_1c_1Vpp_37_20_21_SBM_57_30_12_12_r1` → fields |
| `SBM_ReadBinHeader` | header by key + file-size check (R × C × N_F A-lines) |
| `SBM_FindSidecar` | finds the matching `Info` (next to the `.bin`, by time stamp, LAST, or rebuilt) |
| `SBM_FrameCplx` | complex OCT signal of one frame (math of `MMB_ComputeCplx`) |
| `SBM_SurfaceBand`, `SBM_SurfaceMotion` | Loupas surface velocity (or displacement) across frames |
| `SBM_ToFrames` | record order → [N_B × N_A × N_F] spatial frames (y ascending) |
| `SBM_HarmonicRecovery` | bin N_P phasor, φ_acq compensation, recovered movie, frame spectrum + SNR |
| `SBM_EstimateSource` | needle tip from the cylindrical-wave coherence |
| `SBM_DirectionalFilterOutward2D` | outward-propagation filter of the 2D phasor |
| `SBM_InterpolateSquare` | thesis interpolation to the ≥ 100 × 100 square grid |
| `SBM_MapStats` | median, Ep, Ea (c_ref = 10.54 m/s), valid pixels |
