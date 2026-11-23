# Phase-Gradient Framework for Experimental Multimeridian MB-mode OCE

Author: Estéfano Machado Ferrer

Experimental counterpart of the phase-speed estimation module of the thesis
(Sections 2.3 and 3.6). It runs the four strategies on real data:

| Strategy | Pipeline |
|---|---|
| 1D PG    | phasor at f_exc → local line fit along each meridian → en face interpolation of c |
| 1D PG-TV | same, but interpolate \|k\| → TV (Chambolle–Pock) on \|k\| → c = ω/k_TV |
| 2D PG    | Delaunay interpolation of every sparse en face frame → phasor → local plane fit |
| 2D PG-TV | 2D PG → TV on \|k\| → c = ω/k_TV |

All maps are R²-masked (default R² ≥ 0.7), as in the thesis.

## Quick start

1. Run `OCE_Processing_LiveScript_MBMultimer_v1_0` (or `PP_MBmode.m`) up to
   `MotionAnalysisBorder` → `loaded_phases_Border`.
2. Run `Process_MB_PhaseGradient_RealData.m` section by section.

Or directly:

```matlab
Results = OCE_PG_ProcessMB(loaded_phases_Border, 200e3, 1000, 10, 100, 16);
OCE_PG_PlotResults(Results);
```

## Files

| File | Purpose |
|---|---|
| `OCE_PG_ProcessMB` | Driver: runs all four estimators |
| `OCE_PG1D_AlongMeridians` | 1D PG and 1D PG-TV |
| `OCE_MB_InterpolateEnFaceFrames` | Sparse frames → dense en face movie (2D PG input) |
| `OCE_PG2D` | 2D PG and 2D PG-TV (works on any en face movie, e.g. raster) |
| `OCE_TVDenoise_ChambollePock2D` | TV solver (`'neumann'` default, `'legacy'` = simulation code) |
| `OCE_HarmonicPhasor` | Phasor at f_exc (`'dtft'` exact, `'fftbin'` = simulation code) |
| `OCE_MB_MeridianGeometry` / `OCE_MB_ReshapeToMeridians` / `OCE_MB_EnFaceGrid` | Acquisition geometry |
| `OCE_EstimateDominantFrequency` | Checks that f_exc matches the measured motion |
| `OCE_PG_PlotResults`, `OCE_PG_ROIStats` | Figures, Ep and CNR |
| `Demo_OCE_PG_Synthetic` | Validation with known ground truth |

## Conventions

- Maps are `[Ny x Nx]` (meshgrid, rows = y): `imagesc(R.x_mm, R.y_mm, map); axis xy`.
- Meridians: angles `0:180/N:...`, r from −L/2 to +L/2, every second meridian reversed
  (same as `AcqVectors_MultimeridianMB.m` / `IDXSelection.m`). Verified on the
  25-09-2026 data: the detected surface fits a plane to 0.9 px with this
  convention vs 12–14 px without the reversal.
- The y direction (sign of the meridian angles) could not be verified from the data.
  If your scanner rotates the other way, use `'AngleOffset_deg'`/`'FlipRadial'`.
- `MuTV = []` → `0.1·median|k|`, which matches the thesis choice (μ = 250 rad/m for |k| ≈ 2200–3700 rad/m).

## Differences from the simulation code

- TV: pixels outside the mask are decoupled (Neumann) instead of fixed to 0. The legacy
  behaviour biases |k| toward 0 at mask edges (overestimated speed). Use
  `'TVBoundary','legacy'` to reproduce the thesis numbers.
- Phasor: exact DTFT at f_exc instead of the nearest FFT bin (the records hold only ~4–5
  excitation periods). Use `'PhasorMethod','fftbin'` to reproduce the thesis numbers.
- Interpolation uses the exact (x, y) sample positions, not their nearest reduced-grid pixel.
- Outputs are not cropped; the kernel border is NaN, so all four maps share one grid.
