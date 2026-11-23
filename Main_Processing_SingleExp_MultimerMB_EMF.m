%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% Main_Processing_SingleExp_MultimerMB_EMF
%
% Complete processing of ONE Multimeridian MB-mode OCE acquisition, from
% the raw LabVIEW .bin file up to the four phase-speed elastograms of the
% thesis (1D PG, 1D PG-TV, 2D PG, 2D PG-TV), computed both WITHOUT and
% WITH directional (outward-propagation) filtering.
%
% Based on PP_MBmode.m (same reading, B-mode, surface, Loupas motion and
% masked depth video) and on PP_MBmode_Batch_PG.m (05-10-2026), with the
% 07-10-2026 changes:
%   * narrowband temporal filter: pass band = -3 dB width of the spectral
%     peak closest to f_exc, capped at Cfg.Band.MaxBW_Hz (< 1 kHz), applied
%     with a zero-phase Gaussian filter (MMB_GaussianBandpass);
%   * the phasor and the speed c = 2 pi f/|k| use the measured peak
%     frequency f_an (when the peak is clearly defined);
%   * directional filtering of every meridian (outward waves only,
%     MMB_DirectionalFilterOutward) as an extra branch for comparison;
%   * turbo colormap (thesis) and percentile-based color limits;
%   * particle velocity in mm/s; masked all-meridian depth video.
%
% Run it directly (edit ExpName in section 0) or from
% Main_Processing_Batch_MultimerMB_EMF.m, which sets BatchMode = true,
% ExpName and (optionally) CfgOverride before calling this script.
%
% Outputs (Results/<ExpName>/):
%   Bmode_Border, Spectrum_Band, TimeProfile, SpaceTimeMap_noDF_vs_DF,
%   DirectionalFilter_kf, Elastograms_noDF_vs_DF  (.png + .fig)
%   Video_EnFace_noDF_vs_DF.mp4, Video_Depth_AllMeridians_Masked.mp4,
%   Video_Depth_RepMeridian.mp4 (one representative meridian, for slides)
%   Results_MMB.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Configuration

if ~exist('BatchMode', 'var') || ~BatchMode
    clc; close all;
    BatchMode = false;
    ExpName = 'PC_QH_285us_3c_1Vpp_450_50_100_OCE_8_30_12_12_r1';
end

Cfg = struct();

Cfg.RootDir    = fileparts(mfilename('fullpath'));   % repository root
Cfg.DataDir    = fullfile(Cfg.RootDir, 'Data');
Cfg.ResultsDir = fullfile(Cfg.RootDir, 'Results');
Cfg.CodesDir   = fullfile(Cfg.RootDir, 'Codes');

% OCT system (PP_MBmode.m)
Cfg.OCT.a_scan_rate_kHz   = 200;
Cfg.OCT.axial_pixel_res   = 7.1;     % [um]
Cfg.OCT.IDXrefrac         = 1.39;
Cfg.OCT.center_wavelength = 1.300;   % [um]

% Sizing (PP_MBmode.m)
Cfg.Cut_Depth = [30 750];
Cfg.Cut_Time  = [1 449];
Cfg.i_thresh  = [70 130];
Cfg.PeakThresMult    = 10;           % FindSurface
Cfg.PeakThresWinSize = 30;
Cfg.BmodeMaskThres   = 75;           % [dB]

% Surface motion (PP_MBmode.m: Loupas 10 px, Delta 5, lowess 5 %)
Cfg.Motion.MotionType      = 2;      % particle velocity
Cfg.Motion.Delta           = 5;
Cfg.Motion.LoupasAxialWin  = 10;
Cfg.Motion.SmoothingWinPer = 0.05;

% Narrowband temporal filter
Cfg.Band.SearchRel        = 0.20;    % peak search window f_exc*(1 -/+ 0.2)
Cfg.Band.MaxBW_Hz         = 600;     % cap of the -3 dB width
Cfg.Band.MinProminence_dB = 3;
Cfg.Band.UsePeakFrequency = false;   % false: f_an = f_exc (drive frequency); true: measured peak

% Excitation point (needle). It is located from the data (MMB_EstimateSource)
% because in the 06-10-2026 acquisitions it is ~1-2 mm away from the centre
% of the meridian pattern. Positions closer than NeedleRadius_mm to it
% (needle shaft + near field) are blanked, and the same disk is removed
% from the elastograms (previously: +-9 positions around each meridian
% centre = 1.15 mm).
Cfg.Source.Estimate       = true;
Cfg.Source.MinCoherence   = 0.30;    % below this, fall back to the scan centre
% Session-wide needle position: the needle was not moved during the
% 06-10-2026 session (all files acquired 01:01-01:29). The joint maximum of
% the coherence maps of all >= 4-meridian acquisitions
% (MMB_SourceCoherenceMap, summed) is (0.70, 1.10) mm. Single-acquisition
% fits agree at 2.9-4.5 kHz but wander at 2 kHz (lambda ~ 5.5 mm), so the
% fixed session value is used; the per-experiment fit is still reported.
% Set Fixed_mm = [] to use the per-experiment fit instead.
Cfg.Source.Fixed_mm       = [0.70 1.10];
% !! SESSION-SPECIFIC (06-10-2026, old ~1 mm needle). For a new session /
% needle: set Fixed_mm = [] (or the new session's joint fit) and
% NeedleRadius_mm below. For the 3.54 mm FWHM needle (20-10-2026) the
% simulations show the 2D PG bias is confined to r < ~4.25 mm (2.83 sigma),
% so use NeedleRadius_mm ~ 4.25 there.
Cfg.NeedleRadius_mm       = 1.2;
Cfg.ShaftZ                = 2.5;       % robust z-score of noisy (shaft) columns
Cfg.ShaftSearch_mm        = 3;

% Directional filter (outward propagation along each meridian)
Cfg.DF.SigmaKappaRel = 0.5;
Cfg.DF.SigmaFreqBins = 2;
Cfg.DF.PadFactor     = 4;

% Phase-gradient estimators (thesis: 11 px kernel, 100 x 100 grid, R2 0.7)
Cfg.PG.KernelSize_px = 11;
Cfg.PG.GridSize      = 100;
Cfg.PG.R2Threshold   = 0.70;
Cfg.PG.SpeedLimits   = [0.5 40];
Cfg.PG.MuTVRel       = 0.1;
Cfg.PG.TVIterations  = 500;

% Figures and videos
Cfg.ShowFigures     = ~BatchMode;
Cfg.MakeVideos      = true;
Cfg.MakeDepthVideo  = true;
Cfg.Video.Time_ini  = 10;
Cfg.Video.FrameRate = 20;
Cfg.Video.Quality   = 60;
Cfg.CLimPct         = 99;            % (kept for compatibility)
Cfg.WaveCLimPct     = 97;            % space-time maps and en face videos: +-P97 of |v|, blue-white-red
% Depth video over the B-mode (20-10-2026): blue-white-red, transparent at
% v = 0, limits = P90 of |v| over pixels with B-mode >= 90 dB (the phase of
% low-signal pixels is noise and used to set the scale to ~2.4 mm/s while
% the wave is ~0.15 mm/s rms)
Cfg.DepthVideo.CLimPct     = 90;
Cfg.DepthVideo.MinBmode_dB = 90;

% Overrides from the batch script
if exist('CfgOverride', 'var') && isstruct(CfgOverride)
    fn = fieldnames(CfgOverride);
    for i = 1:numel(fn)
        Cfg.(fn{i}) = CfgOverride.(fn{i});
    end
end

addpath(genpath(Cfg.CodesDir));
warning('off', 'signal:findpeaks:largeMinPeakHeight');
if ~Cfg.ShowFigures
    set(0, 'DefaultFigureVisible', 'off');
end

tStartAll = tic;

%% 1. Experiment parameters (file name + .bin header)

E = MMB_ParseExperimentName(ExpName);

hit = dir(fullfile(Cfg.DataDir, '**', [ExpName '.bin']));
if isempty(hit)
    error('Raw file %s.bin not found under %s', ExpName, Cfg.DataDir);
end
binFile = fullfile(hit(1).folder, hit(1).name);

[ScanInfo, headerBytes, FileCheck] = MMB_ReadBinHeader(binFile);
if ~FileCheck.IsComplete
    error(['%s is incomplete: %.2f of %.2f GB on disk (%d of %d meridians). ' ...
           'The acquisition was interrupted; the file is not processed.'], ...
           ExpName, FileCheck.ActualBytes/1e9, FileCheck.ExpectedBytes/1e9, ...
           FileCheck.CompleteMeridians, ScanInfo.No_3Dscans);
end

nS   = ScanInfo.samples_in_Aline;
nA   = ScanInfo.Alines_in_Bframe;
nPos = ScanInfo.Bframes_in_3Dscan;
nMer = ScanInfo.No_3Dscans;
nTot = nPos*nMer;
scanLength_mm = ScanInfo.Ver_scan_length_mm;     % 12 mm (header now matches the log)

fs = Cfg.OCT.a_scan_rate_kHz*1e3;
f_exc = E.f_exc;

% Particle velocity scale: v = lambda0/(4 pi n) * dphi * fs  [mm/s per rad]
VelScale = Cfg.OCT.center_wavelength*1e-3/(4*pi*Cfg.OCT.IDXrefrac)*fs;

outDir = fullfile(Cfg.ResultsDir, ExpName);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('%s | %d meridians x %d positions, %d A-lines, L = %g mm, f_exc = %.1f Hz\n', ...
    ExpName, nMer, nPos, nA, scanLength_mm, f_exc);

%% 2. Raw data (memory-mapped) and DC spectrum (ReadCplx_Volume_BS)

m = memmapfile(binFile, 'Offset', headerBytes, 'Format', {'uint16', [nS nA nTot], 'raw'});

tIdx = Cfg.Cut_Time(1):min(Cfg.Cut_Time(2), nA);
zIdx = Cfg.Cut_Depth(1):min(Cfg.Cut_Depth(2), nS/2);
nT = numel(tIdx);
nZ = numel(zIdx);
hannK = hann(nS);

Fringes_DC = zeros(nS, nTot);
for pos = 1:nTot
    Fringes_DC(:,pos) = median(double(m.Data.raw(:,:,pos)), 2);
end
DC_spectrum = median(Fringes_DC, 2);
clear Fringes_DC

%% 3. B-mode (one meridian at a time)

Bmode = zeros(nZ, nTot);
for b = 1:nMer
    posIdx = (b-1)*nPos + (1:nPos);
    C = MMB_ComputeCplx(m, posIdx, DC_spectrum, hannK, tIdx, zIdx);
    Bmode(:, posIdx) = mean(abs(C), 3).';
end
Bmode_IntLog = real(20*log10(Bmode));

Xaxis = linspace(0, scanLength_mm*nMer, nTot);
Zaxis = (0:nZ-1)'*Cfg.OCT.axial_pixel_res*1e-3/Cfg.OCT.IDXrefrac;
Time  = (0:nT-2)/fs*1e3;           % [ms] (velocity = diff of phase)

%% 4. Surface detection and B-mode mask (PP_MBmode.m)

OCE_system = struct('NumLateralPos', nTot, 'NumMrept', nA, 'DepthSize', nS/2, ...
    'i_thresh_low', Cfg.i_thresh(1), 'i_thresh_high', Cfg.i_thresh(2), ...
    'Cut_Depth_ini', zIdx(1), 'Cut_Depth_end', zIdx(end), 'NewDepthSize', nZ, ...
    'Cut_Time_ini', tIdx(1), 'Cut_Time_end', tIdx(end), 'NewTimeSize', nT, ...
    'PeakThresMult', Cfg.PeakThresMult, 'PeakThresWinSize', Cfg.PeakThresWinSize, ...
    'Jump_Lat_pos', 0);
OCT_system = struct('NumAngles', nMer, 'a_scan_rate', Cfg.OCT.a_scan_rate_kHz, ...
    'axial_pixel_res', Cfg.OCT.axial_pixel_res, 'IDXrefrac', Cfg.OCT.IDXrefrac, ...
    'x_scan_width', scanLength_mm, 'spec_len', nS);

Border = FindSurface(Bmode, OCE_system, OCT_system);
Border.Idx(end) = nan;
Border.Idx = medfilt2(Border.Idx, [1 7], 'Symmetric');
Border.Idx = smooth(Border.Idx, 0.005, 'lowess');
Border.Idx = Border.Idx - 4;
Border.Idx(Border.Idx < 20) = 20;
Border.Idx(Border.Idx > nZ-1-50) = nZ-1-50;
Border.Idx = Border.Idx + 4;

Bmode_Mask = Bmode_IntLog > Cfg.BmodeMaskThres;

zmm = @(idx) (idx-1)*Cfg.OCT.axial_pixel_res*1e-3/Cfg.OCT.IDXrefrac;
fig = figure('Position', [40 40 1500 560], 'Color', 'w');
imagesc(Xaxis, Zaxis, Bmode_IntLog); clim(Cfg.i_thresh); colormap(gray); hold on;
plot(Xaxis, zmm(Border.Idx), 'r', 'LineWidth', 1);
for b = 1:nMer-1, xline(scanLength_mm*b, 'y:'); end
hold off;
xlabel('x-axis (mm), meridians concatenated'); ylabel('z-axis (mm)');
title(sprintf('B-mode and detected surface - %d meridians, %s', nMer, E.Repetition));
set(gca, 'FontSize', 14);
MMB_SaveFigure(fig, fullfile(outDir, 'Bmode_Border'));

%% 5. Surface particle velocity (MotionAnalysisBorder, Loupas) + depth velocity

Options = Cfg.Motion;
Options.LineOrFrame = 1;

loaded_phases_Border = nan(nTot, nT-1);
VelDepthRaw = cell(1, nMer);

for b = 1:nMer
    posIdx = (b-1)*nPos + (1:nPos);
    C = MMB_ComputeCplx(m, posIdx, DC_spectrum, hannK, tIdx, zIdx);
    [~, V] = evalc('MotionAnalysisBorder(C, round(Border.Idx(posIdx)), Options)');
    loaded_phases_Border(posIdx,:) = V;
    if Cfg.MakeDepthVideo
        % Depth-resolved particle velocity (MotionAnalysis, frame mode, no Loupas)
        VelDepthRaw{b} = single(diff(unwrap(-angle(C), [], 3), 1, 3));   % [pos x z x t]
    end
end
clear C

loaded_phases_Border = loaded_phases_Border*VelScale;      % [mm/s]
validRow = all(isfinite(loaded_phases_Border), 2);

%% 6. Spectrum and narrowband selection (pass 1: all valid traces)

Band = MMB_SelectNarrowBand(loaded_phases_Border(validRow, :), fs, f_exc, ...
    'SearchRel', Cfg.Band.SearchRel, 'MaxBW_Hz', Cfg.Band.MaxBW_Hz, ...
    'MinProminence_dB', Cfg.Band.MinProminence_dB);
if ~Cfg.Band.UsePeakFrequency
    Band.f_an = f_exc;
    Band.f_lo = f_exc - Band.BW_Hz/2;
    Band.f_hi = f_exc + Band.BW_Hz/2;
end
f_an = Band.f_an;

fprintf('    pass 1: peak %.0f Hz (clear: %d, prominence %.1f dB), -3 dB width %.0f Hz -> band %.0f-%.0f Hz, f_an = %.0f Hz\n', ...
    Band.f_peak, Band.PeakIsClear, Band.Prominence_dB, Band.BW3dB_Hz, Band.f_lo, Band.f_hi, f_an);

%% 7. Zero-phase narrowband filtering

Frames = MMB_GaussianBandpass(loaded_phases_Border, fs, f_an, Band.BW_Hz, 2);
Frames(~validRow, :) = NaN;

%% 7b. Excitation point (needle) location and near-field mask
% Pass 1 (above): band from all valid traces. The needle itself gives
% broadband, high-amplitude traces, so once the source is located the band
% is selected again without the near-field traces (pass 2) and the
% surface motion is filtered again.

Geom = OCE_MB_MeridianGeometry(nPos, nMer, scanLength_mm, 'ReverseAlternate', true);

Src = struct('x0_mm', 0, 'y0_mm', 0, 'Offset_mm', 0, 'k_radpermm', NaN, 'c_mps', NaN, ...
    'Coherence', NaN, 'Coherence_Centered', NaN, 'c_Centered_mps', NaN, 'Used', 'centre');
if Cfg.Source.Estimate
    S_ = MMB_EstimateSource(Frames, Geom, fs, f_an);
    if S_.Coherence >= Cfg.Source.MinCoherence
        Src = S_;
        Src.Used = 'estimated';
    else
        Src.Coherence = S_.Coherence;
    end
end
Src.Fit_mm = [Src.x0_mm Src.y0_mm];           % per-experiment fit (reported)
if ~isempty(Cfg.Source.Fixed_mm)
    Src.x0_mm = Cfg.Source.Fixed_mm(1);
    Src.y0_mm = Cfg.Source.Fixed_mm(2);
    Src.Used = 'session';
end
Source_mm = [Src.x0_mm Src.y0_mm];

distSrc = hypot(Geom.X_mm(:) - Source_mm(1), Geom.Y_mm(:) - Source_mm(2));
needleTip = distSrc < Cfg.NeedleRadius_mm;

% Needle shaft: the needle enters obliquely, so each meridian crosses the
% shaft at a different point. Those columns carry large broadband motion:
% robust z-score of log(std of the raw trace) per meridian > ShaftZ, only
% within ShaftSearch_mm of the centre, dilated by +-2 positions.
sdRaw = reshape(std(loaded_phases_Border, 0, 2), nPos, nMer);
Ls = log(sdRaw);
zs = (Ls - median(Ls, 1, 'omitnan'))./(1.4826*mad(Ls, 1, 1));
% The needle body above the phantom is also detected as a "surface" much
% shallower than its neighbours: |Border - median filter(31)| > 30 px.
bRef = medfilt1(Border.Idx(:), 31, 'omitnan', 'truncate');
spike = reshape(abs(Border.Idx(:) - bRef) > 30, nPos, nMer);
shaft = (zs > Cfg.ShaftZ | spike) & abs(Geom.s_mm) < Cfg.ShaftSearch_mm;
shaft = imdilate(shaft, ones(5, 1));
needleShaft = shaft(:);

needle = needleTip | needleShaft;

% Pass 2: narrowband selection and filtering without the near-field traces
Band = MMB_SelectNarrowBand(loaded_phases_Border(validRow & ~needle, :), fs, f_exc, ...
    'SearchRel', Cfg.Band.SearchRel, 'MaxBW_Hz', Cfg.Band.MaxBW_Hz, ...
    'MinProminence_dB', Cfg.Band.MinProminence_dB);
if ~Cfg.Band.UsePeakFrequency
    Band.f_an = f_exc;
    Band.f_lo = f_exc - Band.BW_Hz/2;
    Band.f_hi = f_exc + Band.BW_Hz/2;
end
f_an = Band.f_an;
Frames = MMB_GaussianBandpass(loaded_phases_Border, fs, f_an, Band.BW_Hz, 2);
Frames(~validRow, :) = NaN;

fprintf('    pass 2: peak %.0f Hz (clear: %d), -3 dB width %.0f Hz -> band %.0f-%.0f Hz, f_an = %.0f Hz\n', ...
    Band.f_peak, Band.PeakIsClear, Band.BW3dB_Hz, Band.f_lo, Band.f_hi, f_an);
fprintf('    source at (%.2f, %.2f) mm [%s; own fit (%.2f, %.2f)], cylindrical fit: c = %.2f m/s, coherence %.2f (centred: %.2f)\n', ...
    Src.x0_mm, Src.y0_mm, Src.Used, Src.Fit_mm(1), Src.Fit_mm(2), Src.c_mps, Src.Coherence, Src.Coherence_Centered);

% Spectrum and selected band (pass 2)
sigmaG = (Band.BW_Hz/2)/sqrt(log(2));
Hg = exp(-(Band.freqs - f_an).^2/(2*sigmaG^2));
fig = figure('Position', [40 40 1100 520], 'Color', 'w');
plot(Band.freqs/1e3, Band.P_dB, 'k', 'LineWidth', 1.4); hold on;
plot(Band.freqs/1e3, 20*log10(Hg + eps), 'Color', [0 0.45 0.74], 'LineWidth', 1.4);
yl = [-40 3];
patch([Band.f_lo Band.f_hi Band.f_hi Band.f_lo]/1e3, yl([1 1 2 2]), [0.3 0.75 0.93], ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');
xline(f_exc/1e3, 'r--', 'LineWidth', 1.4);
if isfinite(Band.f3dB_lo)
    xline([Band.f3dB_lo Band.f3dB_hi]/1e3, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
end
yline(-3, ':', 'Color', [0.5 0.5 0.5]);
hold off; grid on; box on;
xlim([0 min(2.2*f_exc, fs/2)]/1e3); ylim(yl);
xlabel('Frequency [kHz]'); ylabel('Mean power [dB, re. peak]');
legend({'Mean spectrum (surface)', 'Gaussian filter |H|', ...
    sprintf('Pass band %.0f-%.0f Hz', Band.f_lo, Band.f_hi), ...
    sprintf('f_{exc} = %.0f Hz', f_exc), sprintf('-3 dB width %.0f Hz', Band.BW3dB_Hz)}, ...
    'Location', 'northeast');
title(sprintf('Narrowband selection - f_{exc} = %.0f Hz, peak = %.0f Hz', f_exc, Band.f_peak));
set(gca, 'FontSize', 13);
MMB_SaveFigure(fig, fullfile(outDir, 'Spectrum_Band'));

% Mask of the depth video: tissue below the detected surface, no needle
% (needle columns dilated by +-3 positions so that the shaft does not set
% the color scale of the video)
needleV = imdilate(reshape(needle, nPos, nMer), ones(7, 1));
needleV = needleV(:);
VideoMask = Bmode_Mask;
for c = 1:nTot
    bI = Border.Idx(c);
    if isfinite(bI) && ~needleV(c)
        VideoMask(1:max(round(bI)-1, 0), c) = false;
    else
        VideoMask(:, c) = false;
    end
end

repPos = round(0.3*nPos);
cand = find(validRow & ~needle);
[~, k] = min(abs(cand - repPos));
repPos = cand(k);

fig = figure('Position', [40 40 1100 480], 'Color', 'w');
plot(Time, loaded_phases_Border(repPos,:), 'Color', [0.2 0.2 0.2], 'LineWidth', 1); hold on;
plot(Time, Frames(repPos,:), 'r', 'LineWidth', 1.6);
hold off; grid on;
xlabel('time (ms)'); ylabel('Particle velocity [mm/s]');
legend({'Raw (Loupas, surface)', sprintf('Narrowband %.0f-%.0f Hz (zero phase)', Band.f_lo, Band.f_hi)});
title(sprintf('Time signal at position %d (meridian %d)', repPos, ceil(repPos/nPos)));
xlim([Time(1) Time(end)]);
set(gca, 'FontSize', 13);
MMB_SaveFigure(fig, fullfile(outDir, 'TimeProfile'));

%% 8. Directional filtering (outward propagation)

FramesIn = Frames;
FramesIn(needle, :) = NaN;

[FramesDF, DFInfo] = MMB_DirectionalFilterOutward(FramesIn, Geom, fs, f_an, ...
    'SigmaKappaRel', Cfg.DF.SigmaKappaRel, 'SigmaFreqBins', Cfg.DF.SigmaFreqBins, ...
    'PadFactor', Cfg.DF.PadFactor, 'Source_mm', Source_mm);

fprintf('    directional filter: k-f peak speed %.2f m/s, outward energy %.0f %% (median over meridians)\n', ...
    DFInfo.c_kf_mps, 100*median(DFInfo.OutwardEnergyFraction));

% k-f spectrum of one meridian before / after
fig = figure('Position', [40 40 1200 480], 'Color', 'w');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
kk = DFInfo.kappa_cyc_per_mm; ff = DFInfo.f_Hz/1e3;
spB = 20*log10(DFInfo.Spec_Before/max(DFInfo.Spec_Before(:)));
spA = 20*log10(DFInfo.Spec_After/max(DFInfo.Spec_Before(:)));
ttls = {'Before (meridian 1)', 'After outward filter'};
spc = {spB, spA};
for i = 1:2
    ax = nexttile(tl);
    imagesc(ax, kk, ff, spc{i}.'); axis(ax, 'xy');
    colormap(ax, turbo); clim(ax, [-40 0]);
    xlim(ax, 4*DFInfo.kappa0_cyc_per_m*1e-3*[-1 1]); ylim(ax, 2*f_an/1e3*[-1 1]);
    xlabel(ax, '\kappa [cycles/mm]'); ylabel(ax, 'f [kHz]'); title(ax, ttls{i});
    set(ax, 'FontSize', 12);
end
cb = colorbar(ax); ylabel(cb, '|V(\kappa,f)| [dB]');
title(tl, sprintf('Wavenumber-frequency spectrum, f_{an} = %.0f Hz', f_an), 'FontWeight', 'bold');
MMB_SaveFigure(fig, fullfile(outDir, 'DirectionalFilter_kf'));

%% 9. Space-time maps (without / with directional filtering)

SpaceTime   = FramesIn   - mean(FramesIn,   2, 'omitnan');
SpaceTimeDF = FramesDF   - mean(FramesDF,   2, 'omitnan');
% Blue-white-red (white = 0), limits +-P97 of |v| of the no-DF map, needle in gray
[fig, CLst] = MMB_PlotSpaceTime({SpaceTime, SpaceTimeDF}, {'Narrowband filtered', 'Narrowband + directional (outward) filter'}, ...
    Xaxis, Time, scanLength_mm, sprintf('Space-time maps - %d meridians, f_{an} = %.0f Hz', nMer, f_an), ...
    'CLimPct', Cfg.WaveCLimPct);
MMB_SaveFigure(fig, fullfile(outDir, 'SpaceTimeMap_noDF_vs_DF'));

% '_Full' version: same data with the needle near field / shaft KEPT (DF
% re-run on the full traces) and the SAME limits, so the near field
% saturates while the far field reads exactly as in the cropped figure
FramesDFFull = MMB_DirectionalFilterOutward(Frames, Geom, fs, f_an, ...
    'SigmaKappaRel', Cfg.DF.SigmaKappaRel, 'SigmaFreqBins', Cfg.DF.SigmaFreqBins, ...
    'PadFactor', Cfg.DF.PadFactor, 'Source_mm', Source_mm);
fig = MMB_PlotSpaceTime({Frames - mean(Frames, 2, 'omitnan'), FramesDFFull - mean(FramesDFFull, 2, 'omitnan')}, ...
    {'Narrowband filtered (near field kept)', 'Narrowband + directional (outward) filter (near field kept)'}, ...
    Xaxis, Time, scanLength_mm, sprintf('Space-time maps - %d meridians, f_{an} = %.0f Hz - full, cropped-version limits', nMer, f_an), ...
    'CLim', CLst);
MMB_SaveFigure(fig, fullfile(outDir, 'SpaceTimeMap_noDF_vs_DF_Full'));

%% 10. Masked depth motion video, all meridians (PP_MBmode.m style)

if Cfg.MakeDepthVideo && Cfg.MakeVideos
    VelDepth = zeros(nZ, nTot, nT-1, 'single');
    for b = 1:nMer
        posIdx = (b-1)*nPos + (1:nPos);
        v = MMB_GaussianBandpass(VelDepthRaw{b}*VelScale, fs, f_an, Band.BW_Hz, 3);  % [pos x z x t]
        VelDepth(:, posIdx, :) = permute(v, [2 1 3]);
        VelDepthRaw{b} = [];
    end
    dvArgs = {'FrameIdx', Cfg.Video.Time_ini:nT-1, 'FrameRate', Cfg.Video.FrameRate, ...
        'MinBmode_dB', Cfg.DepthVideo.MinBmode_dB, 'IThresh', Cfg.i_thresh};
    CLdv = MMB_WriteDepthVideo(VelDepth, Bmode_IntLog, VideoMask, Xaxis, Zaxis, Time, ...
        fullfile(outDir, 'Video_Depth_AllMeridians_Masked.mp4'), dvArgs{:}, ...
        'Quality', 50, 'CLimPct', Cfg.DepthVideo.CLimPct, ...
        'MeridianLength_mm', scanLength_mm, 'MainTitle', sprintf('%d mer, %.0f Hz |', nMer, f_an));

    % '_Full': needle columns kept (tissue below the surface only), same CL
    MaskFull = Bmode_Mask;
    for c = 1:nTot
        if isfinite(Border.Idx(c))
            MaskFull(1:max(round(Border.Idx(c))-1, 0), c) = false;
        else
            MaskFull(:, c) = false;
        end
    end
    MMB_WriteDepthVideo(VelDepth, Bmode_IntLog, MaskFull, Xaxis, Zaxis, Time, ...
        fullfile(outDir, 'Video_Depth_AllMeridians_Masked_Full.mp4'), dvArgs{:}, 'Quality', 50, 'CLim', CLdv, ...
        'MeridianLength_mm', scanLength_mm, 'MainTitle', sprintf('%d mer, %.0f Hz | near field kept |', nMer, f_an));

    % Representative meridian (used in the slides): most reliable pixels,
    % ties -> closest to the needle tip; drawn along s with the tip marked
    [bRep, RepMeridian] = MMB_SelectRepresentativeMeridian(VideoMask, Bmode_IntLog, Geom, Source_mm, ...
        Cfg.DepthVideo.MinBmode_dB);
    posRep = (bRep-1)*nPos + (1:nPos);
    CLrep = MMB_WriteRepMeridianDepthVideo(VelDepth(:, posRep, :), Bmode_IntLog(:, posRep), VideoMask(:, posRep), ...
        Geom, bRep, Source_mm, Zaxis, Time, fullfile(outDir, 'Video_Depth_RepMeridian.mp4'), ...
        dvArgs{:}, 'Quality', 60, 'CLimPct', Cfg.DepthVideo.CLimPct, ...
        'MainTitle', sprintf('%d mer, meridian %d, %.0f Hz |', nMer, bRep, f_an));
    MMB_WriteRepMeridianDepthVideo(VelDepth(:, posRep, :), Bmode_IntLog(:, posRep), MaskFull(:, posRep), ...
        Geom, bRep, Source_mm, Zaxis, Time, fullfile(outDir, 'Video_Depth_RepMeridian_Full.mp4'), ...
        dvArgs{:}, 'Quality', 60, 'CLim', CLrep, ...
        'MainTitle', sprintf('%d mer, meridian %d, %.0f Hz | near field kept |', nMer, bRep, f_an));
    save(fullfile(outDir, 'RepMeridian.mat'), 'RepMeridian');
    clear VelDepth
end
clear VelDepthRaw

%% 11. Phase-gradient estimators: 1D PG, 1D PG-TV, 2D PG, 2D PG-TV

pgArgs = {'KernelSize_px', Cfg.PG.KernelSize_px, 'GridSize', Cfg.PG.GridSize, ...
    'R2Threshold', Cfg.PG.R2Threshold, 'SpeedLimits', Cfg.PG.SpeedLimits, ...
    'ExcludeCenterRadius_mm', 0, ...
    'MuTVRel', Cfg.PG.MuTVRel, 'TVIterations', Cfg.PG.TVIterations, ...
    'KeepMovie', true, 'Verbose', false};

warning('off', 'all');
Res   = OCE_PG_ProcessMB(FramesIn, fs, f_an, scanLength_mm, nPos, nMer, pgArgs{:});
ResDF = OCE_PG_ProcessMB(FramesDF, fs, f_an, scanLength_mm, nPos, nMer, pgArgs{:});
warning('on', 'all');
warning('off', 'signal:findpeaks:largeMinPeakHeight');

% Remove the needle / near-field disk around the excitation point
[XqE, YqE] = meshgrid(Res.x_mm, Res.y_mm);
NearSource = hypot(XqE - Source_mm(1), YqE - Source_mm(2)) < Cfg.NeedleRadius_mm;
names = {'PG1D', 'PG1D_TV', 'PG2D', 'PG2D_TV'};
for i = 1:4
    Res.Maps.(names{i})(NearSource) = NaN;
    ResDF.Maps.(names{i})(NearSource) = NaN;
end

Stats   = MMB_MapStats(Res.Maps);
StatsDF = MMB_MapStats(ResDF.Maps);
disp(table(Stats.Method, Stats.Median, Stats.Ep, Stats.ValidPx, StatsDF.Median, StatsDF.Ep, StatsDF.ValidPx, ...
    'VariableNames', {'Method', 'Median_noDF', 'Ep_noDF', 'Px_noDF', 'Median_DF', 'Ep_DF', 'Px_DF'}));

%% 12. Elastograms (rows: without / with directional filtering)

colT ={'1D PG', '1D PG-TV', '2D PG', '2D PG-TV'};
MapCell = cell(2, 4);
for i = 1:4
    MapCell{1, i} = Res.Maps.(names{i});
    MapCell{2, i} = ResDF.Maps.(names{i});
end
CLimE = MMB_SpeedCLim(MapCell);

fig = MMB_PlotElastogramGrid(MapCell, Res.x_mm, Res.y_mm, {'No DF', 'With DF'}, colT, CLimE, ...
    sprintf('%d Meridians, f_{an} = %.0f Hz, %s', nMer, f_an, E.Repetition), 'RowSuffix', '', 'Source_mm', Source_mm);
MMB_SaveFigure(fig, fullfile(outDir, 'Elastograms_noDF_vs_DF'));

%% 13. En face wavefield video (without / with directional filtering)

if Cfg.MakeVideos
    % Near field of the needle tip removed from the displayed wavefield
    Res.Movie(repmat(NearSource, 1, 1, size(Res.Movie, 3))) = NaN;
    ResDF.Movie(repmat(NearSource, 1, 1, size(ResDF.Movie, 3))) = NaN;
    efArgs = {'FrameIdx', Cfg.Video.Time_ini:size(Res.Movie, 3), 'FrameRate', Cfg.Video.FrameRate, ...
        'Quality', Cfg.Video.Quality, 'MainTitle', sprintf('%d mer, %.0f Hz |', nMer, f_an)};
    CLef = MMB_WriteEnFaceVideo({Res.Movie, ResDF.Movie}, {'No DF', 'With DF (outward)'}, ...
        Res.x_mm, Res.y_mm, Time, Geom, fullfile(outDir, 'Video_EnFace_noDF_vs_DF.mp4'), ...
        'CLimPct', Cfg.WaveCLimPct, efArgs{:});

    % '_Full': movies interpolated from the traces WITH the needle rows,
    % same interpolation as OCE_PG_ProcessMB, same CL (near field saturates)
    MvFull = cell(1, 2); srcFull = {Frames, FramesDFFull};
    for k = 1:2
        MvFull{k} = OCE_MB_InterpolateEnFaceFrames(OCE_MB_ReshapeToMeridians(srcFull{k}, nPos, nMer), Geom, ...
            'GridSize', Cfg.PG.GridSize, 'InterpMethod', 'linear', 'Smoothing', true, 'Sigma_px', 0.8, 'Verbose', false);
    end
    MMB_WriteEnFaceVideo(MvFull, {'No DF (near field kept)', 'With DF (near field kept)'}, ...
        Res.x_mm, Res.y_mm, Time, Geom, fullfile(outDir, 'Video_EnFace_noDF_vs_DF_Full.mp4'), ...
        'CLim', CLef, efArgs{:});
    clear MvFull
end

%% 14. Save

Res   = rmfield(Res, 'Movie');
ResDF = rmfield(ResDF, 'Movie');
DFInfo.Spec_Before = single(DFInfo.Spec_Before);
DFInfo.Spec_After  = single(DFInfo.Spec_After);
DFInfo.Hplus       = single(DFInfo.Hplus);

Info = struct();
Info.Experiment = E;
Info.ScanInfo = ScanInfo;
Info.FileCheck = FileCheck;
Info.ScanLength_mm = scanLength_mm;
Info.fs = fs;
Info.f_exc = f_exc;
Info.f_an = f_an;
Info.Band = rmfield(Band, {'freqs', 'P'});
Info.Spectrum_f = Band.freqs(Band.freqs <= 3*f_exc);
Info.Spectrum_PdB = Band.P_dB(Band.freqs <= 3*f_exc);
Info.VelScale_mmps_per_rad = VelScale;
Info.RepresentativePos = repPos;
Info.Source = Src;
Info.Source_mm = Source_mm;
Info.nNeedleTip = nnz(needleTip);
Info.nNeedleShaft = nnz(needleShaft & ~needleTip);
Info.ValidSurfaceFraction = mean(validRow);
Info.Cfg = Cfg;
Info.ProcessingTime_min = toc(tStartAll)/60;

save(fullfile(outDir, 'Results_MMB.mat'), 'Res', 'ResDF', 'Stats', 'StatsDF', 'Info', 'DFInfo', ...
    'Border', 'needle', 'loaded_phases_Border', 'Frames', 'FramesDF', 'SpaceTime', 'SpaceTimeDF', ...
    'Xaxis', 'Zaxis', 'Time', 'Bmode_IntLog', 'Bmode_Mask', 'VideoMask', 'Geom', '-v7.3');

fprintf('    done in %.1f min | medians noDF: %s | DF: %s m/s\n', Info.ProcessingTime_min, ...
    mat2str(round(Stats.Median', 2)), mat2str(round(StatsDF.Median', 2)));

if ~Cfg.ShowFigures
    close all;
end
