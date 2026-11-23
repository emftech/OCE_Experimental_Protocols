%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% Main_Processing_SingleExp_SBM_EMF
%
% Complete processing of ONE Synchronized BM-mode (SBM) OCE acquisition of
% the Twin OCT, from the raw LabVIEW .bin file (+ its sidecar) up to the
% 2D PG and 2D PG-TV phase-speed elastograms, computed WITHOUT and WITH
% directional (outward-propagation) filtering. It follows section by
% section Main_Processing_SingleExp_MultimerMB_EMF.m (same reading, DC,
% B-mode, FindSurface, Loupas surface motion, figures, videos, Results
% .mat) and the thesis SBM chain (Chapter 5, Sampling_Acq_Sim_SBM_v4):
%
%   MB step                               SBM step (here)
%   M-scan time axis (A-lines)            frame index = equivalent time T_SBM
%   Loupas lag-1 along A-lines            Loupas lag-1 across frames (circular)
%   narrowband Gaussian filter            exact DFT bin N_P (N_F frames = N_P periods)
%   -                                     acquisition-phase compensation (Eq. 59)
%   needle fit on the meridians           coherence fit on the 2D phasor
%   1D outward filter per meridian        2D sector-wise outward filter
%   sparse -> dense interpolation         spline onto >= 100 x 100 square grid
%                                         (Prepare_SBM_2DPhaseGradient)
%   OCE_PG_ProcessMB (1D/2D)              OCE_PG2D (2D PG, 2D PG-TV), as in the thesis
%
% Run it directly (edit ExpName in section 0) or from
% Main_Processing_Batch_SBM_EMF.m (BatchMode = true, ExpName, CfgOverride).
%
% Outputs (SBM_Experiments/Results/<ExpName>/):
%   Bmode_Border, EnFace_Structure, Spectrum_Frames, TimeProfile,
%   Phase_Recovery, DirectionalFilter_k, SpaceTimeMaps_raw_rec_DF,
%   Elastograms_noDF_vs_DF  (.png + .fig)
%   Video_EnFace_noDF_vs_DF.mp4, Video_Depth_CentralBscan.mp4
%   Results_SBM.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Configuration

if ~exist('BatchMode', 'var') || ~BatchMode
    clc; close all;
    BatchMode = false;
    ExpName = 'PC_H_285us_1c_1Vpp_37_20_21_SBM_57_30_12_12_r1';
end

Cfg = struct();

Cfg.RootDir     = fileparts(fileparts(mfilename('fullpath')));   % repository root
Cfg.ExpDir      = fullfile(Cfg.RootDir, 'SBM_Experiments');
Cfg.DataDir     = fullfile(Cfg.ExpDir, 'Data');
Cfg.SidecarDirs = {fullfile(Cfg.ExpDir, 'Sidecars'), Cfg.DataDir};   % SBM_Waveform_*.mat
Cfg.ResultsDir  = fullfile(Cfg.ExpDir, 'Results');
Cfg.CodesDir    = fullfile(Cfg.RootDir, 'Codes');                    % MB codes (reused)
Cfg.SBMCodesDir = fullfile(Cfg.RootDir, 'Codes_SBM_Processing');
Cfg.ProtocolDir = fullfile(Cfg.RootDir, 'SBM_Protocol'); % only to rebuild a missing sidecar

% OCT system (PP_MBmode.m)
Cfg.OCT.a_scan_rate_kHz   = 200;
Cfg.OCT.axial_pixel_res   = 7.1;     % [um]
Cfg.OCT.IDXrefrac         = 1.39;
Cfg.OCT.center_wavelength = 1.300;   % [um]

% Sizing (PP_MBmode.m)
Cfg.Cut_Depth = [30 750];
Cfg.i_thresh  = [70 130];
Cfg.PeakThresMult    = 10;           % FindSurface
Cfg.PeakThresWinSize = 30;
Cfg.BmodeMaskThres   = 75;           % [dB]
Cfg.SurfaceMaxDrop_dB = 20;          % guard: a FindSurface pick more than this below the
                                     % profile maximum (noise peak) is replaced by the first
                                     % peak within this range of the maximum ([] = off)

% Surface motion (PP_MBmode.m: Loupas 10 px, Delta 5) across the frames
Cfg.Motion.MotionType      = 2;      % 2 = particle velocity (as MB), 1 = displacement
Cfg.Motion.Delta           = 5;
Cfg.Motion.LoupasAxialWin  = 10;
Cfg.Motion.SmoothingWinPer = [];     % none: the harmonic is one exact DFT bin
Cfg.Motion.Estimator       = 'band'; % 'band': complex lag-1 correlation over BandDepth_px below
Cfg.Motion.BandDepth_px    = 60;     % the surface (0.3 mm, ~4x less noise on real data);
                                     % 'loupas': MB estimator (LoupasAxialWin, Delta)
Cfg.Motion.BulkCompensation = true;  % remove the bulk (slow vibration) phase between frames,
                                     % moving circular mean over one raster line (P A-lines).
                                     % Essential on real data: 18-10-2026 tests 1.8 -> 0.2 rad

% Temporal band (MB narrowband step, SBM_SelectBand). The piezo has a
% non-flat response (fundamental resonance ~250 Hz), so, as in MB, the band
% is centred on the HIGHEST spectral peak near f_exc (equivalent-time
% spectrum of the frames, Hann, zero-padded), -3 dB width capped at one bin
% (f_exc/N_P = native resolution). 'exact' = bin N_P only (thesis).
Cfg.Band.Mode             = 'peak';
Cfg.Band.SearchRel        = 0.5;     % search window f_exc*(1 -/+ 0.5) (bins N_P/2 .. 3N_P/2)
Cfg.Band.MaxBW_Hz         = [];      % [] = f_exc/N_P (one bin)
Cfg.Band.MinProminence_dB = 3;
Cfg.Band.UsePeakFrequency = true;    % phasor, phi_acq and c = 2 pi f/|k| at the peak frequency
Cfg.Band.Nfft             = 1024;

% Excitation point (needle). Unlike the 06-10 MB session the needle position
% of an SBM session is not known a priori: it is fitted on every acquisition
% (SBM_EstimateSource). Set Fixed_mm = [x0 y0] once a session-wide value is
% known (Main_Processing_Batch_SBM_EMF prints the joint estimate).
Cfg.Source.Estimate     = true;
Cfg.Source.MinCoherence = 0.30;
Cfg.Source.Fixed_mm     = [];
Cfg.NeedleRadius_mm     = 4.25;      % thick needle (20-10-2026, FWHM 3.54 mm): simulated radius beyond
                                     % which 2D PG is within 5 % of c (NearField_Radius.xlsx: 3.75-4.25 mm
                                     % for 220-500 us). Thin needle (18/19-10 data, synthetic test): 1.2
Cfg.FarFieldMargin_mm   = 1;         % far-field ROI r >= NeedleRadius + 1 mm (the first ~1 mm outside
Cfg.RingWidth_mm        = 0.5;       % the disk is biased: interpolation + PG kernel on the edge)
Cfg.ShaftZ              = 2.5;       % robust z-score of noisy (shaft) pixels
Cfg.ShaftSearch_mm      = 3;

% Directional filter (2D, outward from the needle)
Cfg.DF.NSectors      = 8;
Cfg.DF.SigmaKappaRel = 0.5;
Cfg.DF.PadFactor     = 4;

% Phase-gradient estimators (thesis: 11 px kernel, >= 100 x 100 grid, R2 0.7)
Cfg.PG.InterpMethod  = 'spline';     % Prepare_SBM_2DPhaseGradient (thesis runs)
Cfg.PG.GridMin       = 100;
Cfg.PG.KernelSize_px = 11;
Cfg.PG.R2Threshold   = 0.70;
Cfg.PG.SpeedLimits   = [0.5 40];
Cfg.PG.MuTVRel       = 0.1;
Cfg.PG.TVIterations  = 500;

% Figures and videos
Cfg.ShowFigures     = ~BatchMode;
Cfg.MakeVideos      = true;
Cfg.MakeDepthVideo  = true;
Cfg.Video.FrameRate = 10;
Cfg.Video.Quality   = 60;
% Wave-propagation displays (space-time maps, en face and depth videos) use
% the MB convention of 20-10-2026: blue -> white -> red (MMB_RedBlue), zero
% motion white/transparent, NaN mid-gray. Elastograms and k-space spectra
% keep turbo.
Cfg.WaveCLimPct           = 97;      % space-time maps and en face video (P97 of |v|)
Cfg.MakeFullWaveDisplays  = true;    % also "_Full" wave displays WITHOUT the near-field /
                                     % shaft crop, with the SAME colour limits as the cropped
                                     % ones (near field saturates, far field stays readable)
Cfg.DepthVideo.CLimPct    = 90;      % P90 of |v| over reliable pixels (MMB_WriteDepthVideo)
Cfg.DepthVideo.MinBmode_dB = 90;     % overlay only where the B-mode is >= 90 dB (real data
                                     % ~97-115 dB; the synthetic set is ~50-57 dB -> use 40)

if exist('CfgOverride', 'var') && isstruct(CfgOverride)
    fn = fieldnames(CfgOverride);
    for i = 1:numel(fn), Cfg.(fn{i}) = CfgOverride.(fn{i}); end
end

addpath(genpath(Cfg.CodesDir));
addpath(genpath(Cfg.SBMCodesDir));
warning('off', 'signal:findpeaks:largeMinPeakHeight');
if ~Cfg.ShowFigures, set(0, 'DefaultFigureVisible', 'off'); end
tStartAll = tic;

%% 1. Experiment parameters (file name + .bin header + sidecar)

E = SBM_ParseExperimentName(ExpName);
hit = dir(fullfile(Cfg.DataDir, '**', [ExpName '.bin']));
if isempty(hit)
    error('Raw file %s.bin not found under %s', ExpName, Cfg.DataDir);
end
binFile = fullfile(hit(1).folder, hit(1).name);

[ScanInfo, headerBytes, FileCheck] = SBM_ReadBinHeader(binFile);
if ~FileCheck.IsComplete
    error('%s is incomplete: %.2f of %.2f GB (%d of %d frames).', ExpName, ...
        FileCheck.ActualBytes/1e9, FileCheck.ExpectedBytes/1e9, ...
        FileCheck.CompleteFrames, ScanInfo.No_3Dscans);
end
[SInfo, SidecarFile] = SBM_FindSidecar(binFile, E, ScanInfo, Cfg.SidecarDirs, Cfg.ProtocolDir);

nS   = ScanInfo.samples_in_Aline;
R    = SInfo.R;   Crec = SInfo.C;   Nper = R*Crec;
nF   = SInfo.Cfg.nFrames;   nA = SInfo.Cfg.nAscans;   nB = SInfo.Cfg.nBscans;
NP   = SInfo.Cfg.nPeriods;  Nphi = SInfo.Nphi;
T_SBM = SInfo.Ts_SBM;                       % equivalent sampling period [s]
fsEq  = 1/T_SBM;                            % equivalent sampling rate
f_exc = E.f_exc;                            % f_exc = DFT bin N_P exactly (f_an: section 6)
xW = ScanInfo.Hor_scan_length_mm;  yW = ScanInfo.Ver_scan_length_mm;
lam_n = Cfg.OCT.center_wavelength*1e-3/(4*pi*Cfg.OCT.IDXrefrac);   % mm per rad
if Cfg.Motion.MotionType == 2
    VelScale = lam_n/T_SBM;   motionUnits = 'Particle velocity [mm/s]';   motionKind = 'velocity';
else
    VelScale = lam_n*1e6;     motionUnits = 'Displacement [nm]';         motionKind = 'displacement';
end

outDir = fullfile(Cfg.ResultsDir, ExpName);
if ~exist(outDir, 'dir'), mkdir(outDir); end
fprintf('%s | %d x %d A-lines, %d frames (N_P = %d, T_SBM = %.0f us), R = %d x C = %d, f_exc = %.1f Hz\n', ...
    ExpName, nA, nB, nF, NP, T_SBM*1e6, R, Crec, f_exc);
fprintf('    sidecar: %s\n', SidecarFile);

%% 2. Raw data (memory-mapped) and DC spectrum

m = memmapfile(binFile, 'Offset', headerBytes, 'Format', {'uint16', [nS Nper nF], 'raw'});
zIdx = Cfg.Cut_Depth(1):min(Cfg.Cut_Depth(2), nS/2);
nZ = numel(zIdx);
hannK = hann(nS);
DCf = zeros(nS, nF);
for f = 1:nF
    DCf(:, f) = median(double(m.Data.raw(:, :, f)), 2);
end
DC_spectrum = median(DCf, 2);
clear DCf

%% 3. B-mode (frame-averaged, record order) and spatial views

Bmode = zeros(nZ, Nper);
rowC = ceil(nB/2);                                   % central B-scan (for the depth video)
PIc  = flipud(SInfo.PixelIdx);  colsC = PIc(rowC, :); % its A-lines, x ascending
CplxC = complex(zeros(nZ, nA, nF, 'single'));
for f = 1:nF
    C = SBM_FrameCplx(m, f, DC_spectrum, hannK, zIdx);
    Bmode = Bmode + abs(C);
    if Cfg.MakeDepthVideo, CplxC(:, :, f) = single(C(:, colsC)); end
end
Bmode = Bmode/nF;
Bmode_IntLog = real(20*log10(Bmode));
Zaxis = (0:nZ-1)'*Cfg.OCT.axial_pixel_res*1e-3/Cfg.OCT.IDXrefrac;

%% 4. Surface detection and B-mode mask (PP_MBmode.m, along the raster path)

OCE_system = struct('NumLateralPos', Nper, 'NumMrept', nF, 'DepthSize', nS/2, ...
    'i_thresh_low', Cfg.i_thresh(1), 'i_thresh_high', Cfg.i_thresh(2), ...
    'Cut_Depth_ini', zIdx(1), 'Cut_Depth_end', zIdx(end), 'NewDepthSize', nZ, ...
    'PeakThresMult', Cfg.PeakThresMult, 'PeakThresWinSize', Cfg.PeakThresWinSize, 'Jump_Lat_pos', 0);
OCT_system = struct('NumAngles', 1, 'a_scan_rate', Cfg.OCT.a_scan_rate_kHz, ...
    'axial_pixel_res', Cfg.OCT.axial_pixel_res, 'IDXrefrac', Cfg.OCT.IDXrefrac, ...
    'x_scan_width', xW, 'spec_len', nS);
Border = FindSurface(Bmode, OCE_system, OCT_system);
nFixed = 0;
if ~isempty(Cfg.SurfaceMaxDrop_dB)
    B1 = 20*log10(medfilt2(Bmode, [3 3], 'symmetric'));      % same profile as FindSurface
    [pkMax, ~] = max(B1(10:end, :), [], 1);
    for j = 1:Nper
        bj = Border.Idx(j);
        if isfinite(bj) && B1(bj, j) >= pkMax(j) - Cfg.SurfaceMaxDrop_dB, continue; end
        [~, lc] = findpeaks(B1(10:end, j), 'MinPeakHeight', pkMax(j) - Cfg.SurfaceMaxDrop_dB, 'NPeaks', 1);
        if ~isempty(lc), Border.Idx(j) = lc(1) + 9; nFixed = nFixed + 1; end
    end
    clear B1
end
fprintf('    surface: FindSurface, %d of %d A-lines corrected (noise picks > %g dB below the profile peak)\n', ...
    nFixed, Nper, Cfg.SurfaceMaxDrop_dB);
Border.Idx = medfilt2(Border.Idx, [1 7], 'Symmetric');
Border.Idx = smooth(Border.Idx, 0.005, 'lowess')';
Border.Idx = Border.Idx - 4;
Border.Idx(Border.Idx < 20) = 20;
Border.Idx(Border.Idx > nZ-1-50) = nZ-1-50;
Border.Idx = Border.Idx + 4;
Bmode_Mask = Bmode_IntLog > Cfg.BmodeMaskThres;

[BmodeVol, x_mm, y_mm] = SBM_ToFrames(Bmode_IntLog, SInfo, xW, yW);    % [nZ x nB x nA]
SurfDepth = squeeze(SBM_ToFrames(Border.Idx, SInfo, xW, yW));         % [nB x nA]
zmm = @(idx) (idx-1)*Cfg.OCT.axial_pixel_res*1e-3/Cfg.OCT.IDXrefrac;

isAl = false(1, Nper);  isAl(SInfo.PixelIdx(:)) = true;
nShow = min(Nper, 3*(nA + SInfo.Cfg.nTurnSamples));
fig = figure('Position', [40 40 1500 560], 'Color', 'w');
imagesc(1:nShow, Zaxis, Bmode_IntLog(:, 1:nShow)); clim(Cfg.i_thresh); colormap(gray); hold on;
plot(1:nShow, zmm(Border.Idx(1:nShow)), 'r', 'LineWidth', 1);
plot(find(isAl(1:nShow)), zmm(Border.Idx(isAl(1:nShow))), 'g.', 'MarkerSize', 6); hold off;
xlabel('A-line index in the frame record (first 3 B-scan steps; green = acquired pixels)'); ylabel('z-axis (mm)');
title(sprintf('Frame-averaged B-mode and detected surface - %s', E.Repetition)); set(gca, 'FontSize', 14);
MMB_SaveFigure(fig, fullfile(outDir, 'Bmode_Border'));

fig = figure('Position', [40 40 1500 460], 'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');
nexttile(tl); imagesc(x_mm, Zaxis, squeeze(BmodeVol(:, rowC, :))); clim(Cfg.i_thresh); colormap(gca, gray);
hold on; plot(x_mm, zmm(SurfDepth(rowC, :)), 'r.'); hold off;
xlabel('x (mm)'); ylabel('z (mm)'); title(sprintf('Central B-scan (y = %.2f mm)', y_mm(rowC)));
nexttile(tl); imagesc(x_mm, y_mm, squeeze(max(BmodeVol, [], 1))); axis image xy; colormap(gca, gray); colorbar;
xlabel('x (mm)'); ylabel('y (mm)'); title('En face max projection [dB]');
nexttile(tl); imagesc(x_mm, y_mm, zmm(SurfDepth)); axis image xy; colormap(gca, parula); colorbar;
xlabel('x (mm)'); ylabel('y (mm)'); title('Surface depth (mm)');
MMB_SaveFigure(fig, fullfile(outDir, 'EnFace_Structure'));

%% 5. Surface motion across the frames (Loupas, MotionAnalysisBorder style)

Band = [];
for f = 1:nF
    C = SBM_FrameCplx(m, f, DC_spectrum, hannK, zIdx);
    [Band, validBand] = SBM_SurfaceBand(C, Border.Idx, Cfg.Motion, Band, f);
end
clear C
SurfMotion = SBM_SurfaceMotion(Band, validBand, Cfg.Motion) * VelScale;   % [Nper x nF]
clear Band
if Cfg.Motion.BulkCompensation && Cfg.Motion.MotionType == 2
    sdBefore = median(std(SurfMotion, 0, 2, 'omitnan'), 'omitnan') / VelScale;
    [SurfMotion, BulkPhase] = SBM_BulkCompensation(SurfMotion, SInfo.P, VelScale);
    fprintf('    bulk compensation: lag-1 phase SD %.2f -> %.2f rad (median over A-lines)\n', sdBefore, ...
        median(std(SurfMotion, 0, 2, 'omitnan'), 'omitnan') / VelScale);
else
    BulkPhase = [];
end
validRow = all(isfinite(SurfMotion), 2);

[FramesRaw, ~, ~] = SBM_ToFrames(SurfMotion, SInfo, xW, yW);   % [nB x nA x nF], frame order
validPix = all(isfinite(FramesRaw), 3);

%% 6. Harmonic recovery (exact DFT bin N_P + Eq. 59) and frame spectrum

H = SBM_HarmonicRecovery(FramesRaw, SInfo, motionKind, Cfg.Band);
f_an = H.f_an;  binAn = H.BinAn;
if H.Band.Filtered
    fprintf(['    band (%s): peak %.0f Hz (bin %.2f, prominence %.1f dB, -3 dB width %.0f Hz) -> ', ...
        'Gaussian band-pass %.0f +- %.0f Hz, f_an = %.0f Hz%s\n'], Cfg.Band.Mode, H.Band.f_peak, ...
        H.Band.f_peak/(f_exc/NP), H.Band.Prominence_dB, H.Band.BW3dB_Hz, f_an, H.Band.BW_Hz/2, f_an, ...
        ternary(H.Band.PeakIsClear, '', ' (no clear peak: centred at f_exc)'));
end
fprintf('    excitation bin %d: SNR %.1f dB over the other bins (largest non-DC peak at bin %d)\n', ...
    NP, H.SNR_dB, H.PeakBin);
fprintf('    detection: spatial coherence of the compensated phasor at %.0f Hz %.2f vs null P95 %.2f -> %s\n', ...
    f_an, H.Detect.CohExc, H.Detect.Threshold, ternary(H.Detect.Detected, 'EXCITATION DETECTED', 'NOT DETECTED'));
if ~H.Detect.Detected
    warning('SBM:NoExcitation', ['%s: no excitation-locked wave detected at %.0f Hz. ', ...
        'The elastograms below are computed on noise and must not be interpreted.'], ExpName, f_an);
end

fig = figure('Position', [40 40 1100 480], 'Color', 'w');
stem(H.f_equiv_Hz/1e3, 20*log10(H.Spectrum/max(H.Spectrum(2:end))), 'k', 'filled', 'MarkerSize', 3); hold on;
stem(H.f_equiv_Hz(NP+1)/1e3, 0, 'r', 'filled', 'LineWidth', 1.5);
if H.Band.Filtered
    plot(H.Band.freqs/1e3, 10*log10(H.Band.P/max(H.Band.P)), 'b-', 'LineWidth', 1);
    xline([H.Band.f_lo H.Band.f_hi]/1e3, 'b--');  xline(f_an/1e3, 'b-', 'f_{an}', 'LineWidth', 1.5);
    legend({'|DFT| per bin', 'f_{exc} (bin N_P)', 'Hann, zero-padded', 'pass band'}, 'Location', 'northeast');
end
hold off; grid on;
xlabel('Equivalent frequency [kHz] (bin k = k f_{exc}/N_P)'); ylabel('Mean |DFT| [dB re. peak]');
title(sprintf('Frame spectrum: f_{exc} = %.0f Hz on bin N_P = %d, SNR %.1f dB', f_exc, NP, H.SNR_dB));
ylim([-60 3]); set(gca, 'FontSize', 13);
MMB_SaveFigure(fig, fullfile(outDir, 'Spectrum_Frames'));

%% 7. Excitation point (needle) and near-field mask

Src = struct('x0_mm', 0, 'y0_mm', 0, 'Offset_mm', 0, 'k_radpermm', NaN, 'c_mps', NaN, ...
    'Coherence', NaN, 'Coherence_Centered', NaN, 'c_Centered_mps', NaN, 'Used', 'centre');
CohMap = [];
if Cfg.Source.Estimate
    [S_, CohMap] = SBM_EstimateSource(H.Phasor, x_mm, y_mm, f_an);
    if S_.Coherence >= Cfg.Source.MinCoherence
        Src = S_;  Src.Used = 'estimated';
    else
        Src.Coherence = S_.Coherence;
    end
end
Src.Fit_mm = [Src.x0_mm Src.y0_mm];
if ~isempty(Cfg.Source.Fixed_mm)
    Src.x0_mm = Cfg.Source.Fixed_mm(1);  Src.y0_mm = Cfg.Source.Fixed_mm(2);  Src.Used = 'session';
end
Source_mm = [Src.x0_mm Src.y0_mm];

[Xg, Yg] = meshgrid(x_mm, y_mm);
needleTip = SBM_NeedleMask(x_mm, y_mm, Source_mm, Cfg.NeedleRadius_mm);
sdRaw = std(FramesRaw, 0, 3);
Ls = log(sdRaw);
zsc = (Ls - median(Ls(validPix), 'omitnan'))/(1.4826*mad(Ls(validPix), 1));
spike = abs(SurfDepth - medfilt2(SurfDepth, [3 5], 'symmetric')) > 30;     % needle body above the phantom
nearSrc = hypot(Xg - Source_mm(1), Yg - Source_mm(2)) < Cfg.ShaftSearch_mm;
needleShaft = imdilate((zsc > Cfg.ShaftZ | spike) & nearSrc, ones(3));
needle = needleTip | needleShaft | ~validPix;

fprintf('    source at (%.2f, %.2f) mm [%s], cylindrical fit c = %.2f m/s, coherence %.2f (centred %.2f); %d needle/invalid px\n', ...
    Source_mm(1), Source_mm(2), Src.Used, Src.c_mps, Src.Coherence, Src.Coherence_Centered, nnz(needle));

U = H.Phasor;  U(needle) = NaN;
U_full = H.Phasor;  U_full(~validPix) = NaN;       % near field kept (for the "_Full" displays only)

% representative pixel time profile (raw frames vs recovered)
dRep = abs(hypot(Xg - Source_mm(1), Yg - Source_mm(2)) - 3) + 1e3*needle;   % ~3 mm from the needle
[~, iRep] = min(dRep(:));
[bRep, aRep] = ind2sub([nB nA], iRep);
fig = figure('Position', [40 40 1100 480], 'Color', 'w');
plot(H.Time_equiv_ms, squeeze(FramesRaw(bRep, aRep, :)), 'o-', 'Color', [0.3 0.3 0.3]); hold on;
plot(H.Time_equiv_ms, squeeze(H.Recovered(bRep, aRep, :)), 'r-', 'LineWidth', 1.6); hold off; grid on;
xlabel('Equivalent time (m-1) T_{SBM} [ms]'); ylabel(motionUnits);
legend({'Measured (frame order)', 'Harmonic after compensation (Eq. 59)'});
title(sprintf('Pixel (x, y) = (%.2f, %.2f) mm, \\phi_{acq} = %.2f rad', x_mm(aRep), y_mm(bRep), H.PhiAcq(bRep, aRep)));
set(gca, 'FontSize', 13);
MMB_SaveFigure(fig, fullfile(outDir, 'TimeProfile'));

%% 8. Directional filtering (outward from the needle)

[UDF, DFInfo] = SBM_DirectionalFilterOutward2D(U, x_mm, y_mm, Source_mm, ...
    'NSectors', Cfg.DF.NSectors, 'SigmaKappaRel', Cfg.DF.SigmaKappaRel, 'PadFactor', Cfg.DF.PadFactor);
fprintf('    directional filter: kappa0 = %.3f cyc/mm (c = %.2f m/s), outward energy %.0f %%\n', ...
    DFInfo.kappa0_cyc_per_mm, f_an/(DFInfo.kappa0_cyc_per_mm*1e3), 100*DFInfo.OutwardEnergyFraction);
if Cfg.MakeFullWaveDisplays
    UDF_full = SBM_DirectionalFilterOutward2D(U_full, x_mm, y_mm, Source_mm, ...
        'NSectors', Cfg.DF.NSectors, 'SigmaKappaRel', Cfg.DF.SigmaKappaRel, 'PadFactor', Cfg.DF.PadFactor);
end

fig = figure('Position', [40 40 1500 900], 'Color', 'w');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact');
ims = {angle(H.PhasorRaw), H.PhiAcq, angle(U), angle(UDF)};
ttl = {'Sampled (raw) phase', '\phi_{acq} (exact, from the timeline)', 'Recovered phase (Eq. 59)', 'Recovered + outward DF'};
for i = 1:4
    nexttile(tl); h_ = imagesc(x_mm, y_mm, angle(exp(1j*ims{i}))); axis image xy;
    set(h_, 'AlphaData', isfinite(ims{i})); colormap(gca, SBM_PhaseCmap); clim([-pi pi]);
    hold on; plot(Source_mm(1), Source_mm(2), 'w+', 'MarkerSize', 10, 'LineWidth', 2); hold off;
    title(ttl{i}); xlabel('x (mm)'); ylabel('y (mm)');
end
kk = 3*DFInfo.kappa0_cyc_per_mm;
nexttile(tl); imagesc(DFInfo.kx_cyc_per_mm, DFInfo.ky_cyc_per_mm, 20*log10(DFInfo.Spec_Before/max(DFInfo.Spec_Before(:))));
axis image xy; xlim(kk*[-1 1]); ylim(kk*[-1 1]); clim([-40 0]); colormap(gca, turbo);
title('|FFT_2\{U\}| before'); xlabel('\kappa_x [cyc/mm]'); ylabel('\kappa_y [cyc/mm]');
nexttile(tl); imagesc(DFInfo.kx_cyc_per_mm, DFInfo.ky_cyc_per_mm, 20*log10(DFInfo.Spec_After/max(DFInfo.Spec_Before(:))));
axis image xy; xlim(kk*[-1 1]); ylim(kk*[-1 1]); clim([-40 0]); colormap(gca, turbo);
title('after outward DF'); xlabel('\kappa_x [cyc/mm]'); ylabel('\kappa_y [cyc/mm]');
cb = colorbar; cb.Label.String = 'dB';
title(tl, sprintf('Phase compensation and directional filter - f_{exc} = %.0f Hz', f_an), 'FontWeight', 'bold');
MMB_SaveFigure(fig, fullfile(outDir, 'Phase_Recovery'));

%% 9. Interpolation to the square grid and space-time maps

[Uq, xq, yq] = SBM_InterpolateSquare(U, x_mm, y_mm, Cfg.PG.InterpMethod, Cfg.PG.GridMin);
[UDFq, ~, ~] = SBM_InterpolateSquare(UDF, x_mm, y_mm, Cfg.PG.InterpMethod, Cfg.PG.GridMin);
mq = reshape(0:nF-1, 1, 1, []);
Movie   = real(Uq   .* exp(1j*2*pi*binAn*mq/nF));   % recovered movie, equivalent time
MovieDF = real(UDFq .* exp(1j*2*pi*binAn*mq/nF));
Time = H.Time_equiv_ms;

[~, iyc] = min(abs(yq - Source_mm(2)));
[~, ibS] = min(abs(y_mm - Source_mm(2)));
RawRow = squeeze(FramesRaw(ibS, :, :));
STs = {RawRow, squeeze(Movie(iyc, :, :)), squeeze(MovieDF(iyc, :, :))};
xs  = {x_mm, xq, xq};
CLimST = MMB_SymmetricCLim(STs{2}, Cfg.WaveCLimPct);
stTitle = sprintf('Space-time maps - %d x %d SBM, f_{exc} = %.0f Hz', nA, nB, f_an);
fig = plotSpaceTime(STs, xs, Time, CLimST, Source_mm(1), yq(iyc), motionUnits, stTitle, '');
MMB_SaveFigure(fig, fullfile(outDir, 'SpaceTimeMaps_raw_rec_DF'));
if Cfg.MakeFullWaveDisplays                     % same maps without the near-field crop, same CLim
    [Uq_full, ~, ~]   = SBM_InterpolateSquare(U_full,   x_mm, y_mm, Cfg.PG.InterpMethod, Cfg.PG.GridMin);
    [UDFq_full, ~, ~] = SBM_InterpolateSquare(UDF_full, x_mm, y_mm, Cfg.PG.InterpMethod, Cfg.PG.GridMin);
    Movie_full   = real(Uq_full   .* exp(1j*2*pi*binAn*mq/nF));
    MovieDF_full = real(UDFq_full .* exp(1j*2*pi*binAn*mq/nF));
    STf = {RawRow, squeeze(Movie_full(iyc, :, :)), squeeze(MovieDF_full(iyc, :, :))};
    fig = plotSpaceTime(STf, xs, Time, CLimST, Source_mm(1), yq(iyc), motionUnits, stTitle, ...
        ' (near field kept, saturated at the cropped colour limits)');
    MMB_SaveFigure(fig, fullfile(outDir, 'SpaceTimeMaps_raw_rec_DF_Full'));
end

%% 10. Depth motion video of the central B-scan (PP_MBmode.m style)

if Cfg.MakeDepthVideo && Cfg.MakeVideos
    % depth-resolved lag-1 phase across frames (circular), then the same
    % harmonic recovery as the surface (per depth, phi_acq of that row)
    dphi = -angle(CplxC(:, :, [2:nF 1]) .* conj(CplxC));           % [nZ x nA x nF], Loupas sign
    if Cfg.Motion.MotionType == 1
        dphi = -unwrap(angle(CplxC .* conj(mean(CplxC, 3))), [], 3);
    end
    if ~isempty(BulkPhase)                                         % same bulk as the surface
        pix = flipud(SInfo.PixelIdx);
        dphi = angle(exp(1j*(double(dphi) - reshape(BulkPhase(pix(rowC, :), :), 1, nA, nF))));
    end
    ph = flipud(SInfo.PhiAcq);  phiRow = ph(rowC, :) * binAn/NP;
    if Cfg.Motion.MotionType == 2, phiRow = phiRow + pi*binAn/nF; end
    D = double(dphi) - mean(double(dphi), 3);
    Ud = (2/nF)*sum(D .* exp(-1j*2*pi*binAn*mq/nF), 3) .* exp(-1j*phiRow);
    VelDepth = single(real(Ud .* exp(1j*2*pi*binAn*mq/nF)) * VelScale);   % [nZ x nA x nF]
    BmodeC = squeeze(BmodeVol(:, rowC, :));
    MaskC = BmodeC > Cfg.BmodeMaskThres;
    for a = 1:nA
        bI = SurfDepth(rowC, a);
        if isfinite(bI) && ~needle(rowC, a), MaskC(1:max(round(bI)-1, 0), a) = false;
        else, MaskC(:, a) = false; end
    end
    % colour limit of the cropped video, computed as MMB_WriteDepthVideo does
    % (P of |v| over reliable pixels of every 5th frame), reused for "_Full"
    Mrel = MaskC & medfilt2(BmodeC, [5 5], 'symmetric') >= Cfg.DepthVideo.MinBmode_dB;
    sub = 1:5:nF;  smp = VelDepth(:, :, sub);
    CLdepth = MMB_SymmetricCLim(smp(repmat(Mrel, 1, 1, numel(sub))), Cfg.DepthVideo.CLimPct);
    dvArgs = {'FrameIdx', 1:nF, 'FrameRate', Cfg.Video.FrameRate, 'Quality', 50, 'CLim', CLdepth, ...
        'MinBmode_dB', Cfg.DepthVideo.MinBmode_dB, 'IThresh', Cfg.i_thresh, ...
        'XLabel', sprintf('x (mm), central B-scan (y = %.2f mm)', y_mm(rowC)), ...
        'VLines', Source_mm(1), 'VLineLabel', 'needle tip', 'FigWidth', 900};
    MMB_WriteDepthVideo(VelDepth, BmodeC, MaskC, x_mm, Zaxis, Time, ...
        fullfile(outDir, 'Video_Depth_CentralBscan.mp4'), dvArgs{:}, ...
        'MainTitle', sprintf('Central B-scan, %.0f Hz |', f_an));
    if Cfg.MakeFullWaveDisplays                 % needle columns kept (only above-surface crop)
        MaskCf = BmodeC > Cfg.BmodeMaskThres;
        for a = 1:nA
            bI = SurfDepth(rowC, a);
            if isfinite(bI), MaskCf(1:max(round(bI)-1, 0), a) = false; else, MaskCf(:, a) = false; end
        end
        MMB_WriteDepthVideo(VelDepth, BmodeC, MaskCf, x_mm, Zaxis, Time, ...
            fullfile(outDir, 'Video_Depth_CentralBscan_Full.mp4'), dvArgs{:}, ...
            'MainTitle', sprintf('Central B-scan, %.0f Hz, near field kept |', f_an));
    end
    clear VelDepth dphi D Ud
end
clear CplxC

%% 11. Phase-gradient estimators: 2D PG, 2D PG-TV (without / with DF)

pgArgs = {'InputIsPhasor', true, 'KernelSize_px', Cfg.PG.KernelSize_px, ...
    'R2Threshold', Cfg.PG.R2Threshold, 'SpeedLimits', Cfg.PG.SpeedLimits, ...
    'MuTVRel', Cfg.PG.MuTVRel, 'TVIterations', Cfg.PG.TVIterations, 'Verbose', false};
warning('off', 'all');
PG   = OCE_PG2D(Uq,   xq, yq, fsEq, f_an, pgArgs{:});
PGDF = OCE_PG2D(UDFq, xq, yq, fsEq, f_an, pgArgs{:});
warning('on', 'all');  warning('off', 'signal:findpeaks:largeMinPeakHeight');

NearSource = SBM_NeedleMask(xq, yq, Source_mm, Cfg.NeedleRadius_mm);
Res   = struct('x_mm', xq, 'y_mm', yq, 'Maps', struct('PG2D', PG.Elastogram_2DPG_Masked, 'PG2D_TV', PG.Elastogram_2DPG_TV_Masked), ...
               'R2', PG.R2_map, 'K', PG.K_map, 'K_TV', PG.K_TV_map, 'Phasor', Uq);
ResDF = struct('x_mm', xq, 'y_mm', yq, 'Maps', struct('PG2D', PGDF.Elastogram_2DPG_Masked, 'PG2D_TV', PGDF.Elastogram_2DPG_TV_Masked), ...
               'R2', PGDF.R2_map, 'K', PGDF.K_map, 'K_TV', PGDF.K_TV_map, 'Phasor', UDFq);
for nm = {'PG2D', 'PG2D_TV'}
    Res.Maps.(nm{1})(NearSource) = NaN;  ResDF.Maps.(nm{1})(NearSource) = NaN;
end
c_ref = 0.9553*sqrt(0.365e6/(3*1000));        % 10.54 m/s, quasi-static c_R of the PC phantom
Stats   = SBM_MapStats(Res.Maps, c_ref);
StatsDF = SBM_MapStats(ResDF.Maps, c_ref);
disp(table(Stats.Method, Stats.Median, Stats.Ep, Stats.ValidPx, StatsDF.Median, StatsDF.Ep, StatsDF.ValidPx, ...
    'VariableNames', {'Method', 'Median_noDF', 'Ep_noDF', 'Px_noDF', 'Median_DF', 'Ep_DF', 'Px_DF'}));
% far-field ROI (r >= R + margin) and median vs ring radius
ffArgs = {Source_mm, Cfg.NeedleRadius_mm, Cfg.FarFieldMargin_mm, Cfg.RingWidth_mm, c_ref};
[StatsFar,   Ring]   = SBM_FarFieldStats(Res.Maps,   xq, yq, ffArgs{:});
[StatsFarDF, RingDF] = SBM_FarFieldStats(ResDF.Maps, xq, yq, ffArgs{:});
fprintf('    far field (r >= %.2f mm): medians noDF %s | DF %s m/s, Ep noDF %s | DF %s %%\n', ...
    Cfg.NeedleRadius_mm + Cfg.FarFieldMargin_mm, mat2str(round(StatsFar.Median', 2)), ...
    mat2str(round(StatsFarDF.Median', 2)), mat2str(round(StatsFar.Ep', 1)), mat2str(round(StatsFarDF.Ep', 1)));
fig = figure('Position', [40 40 900 480], 'Color', 'w');
rc = (Ring.r_in_mm + Ring.r_out_mm)/2;
plot(rc, Ring.PG2D_Median, 'o-', rc, Ring.PG2D_TV_Median, 's-', rc, RingDF.PG2D_Median, 'o--', ...
    rc, RingDF.PG2D_TV_Median, 's--', 'LineWidth', 1.4); hold on;
yline(c_ref, 'k:', 'c_{ref}');  xline(Cfg.NeedleRadius_mm + Cfg.FarFieldMargin_mm, 'k--', 'far-field ROI');  hold off;
colororder(lines(4)); grid on; set(gca, 'FontSize', 13);
xlabel(sprintf('distance from the needle tip (mm), %.1f mm rings', Cfg.RingWidth_mm)); ylabel('median phase speed (m/s)');
legend({'2D PG', '2D PG-TV', '2D PG + DF', '2D PG-TV + DF'}, 'Location', 'best');
title(sprintf('Median vs ring radius (blanked r < %.2f mm)', Cfg.NeedleRadius_mm));
MMB_SaveFigure(fig, fullfile(outDir, 'RingProfile'));

%% 12. Elastograms (rows: without / with directional filtering)

MapCell = {Res.Maps.PG2D, Res.Maps.PG2D_TV; ResDF.Maps.PG2D, ResDF.Maps.PG2D_TV};
CLimE = MMB_SpeedCLim(MapCell);
fig = MMB_PlotElastogramGrid(MapCell, xq, yq, {'No DF', 'With DF'}, {'2D PG', '2D PG-TV'}, CLimE, ...
    sprintf('SBM %d x %d, %d frames, f_{exc} = %.0f Hz, %s', nA, nB, nF, f_an, E.Repetition), ...
    'RowSuffix', '', 'Source_mm', Source_mm);
MMB_SaveFigure(fig, fullfile(outDir, 'Elastograms_noDF_vs_DF'));

%% 13. En face wavefield video (without / with directional filtering)

if Cfg.MakeVideos
    Movie(repmat(NearSource, 1, 1, nF)) = NaN;  MovieDF(repmat(NearSource, 1, 1, nF)) = NaN;
    MMB_WriteEnFaceVideo({Movie, MovieDF}, {'Recovered', 'Recovered + outward DF'}, xq, yq, Time, [], ...
        fullfile(outDir, 'Video_EnFace_noDF_vs_DF.mp4'), 'FrameIdx', 1:nF, ...
        'FrameRate', Cfg.Video.FrameRate, 'Quality', Cfg.Video.Quality, 'CLimPct', Cfg.WaveCLimPct, ...
        'Units', motionUnits, 'MainTitle', sprintf('SBM %dx%d, %.0f Hz |', nA, nB, f_an));
    if Cfg.MakeFullWaveDisplays
        % same colour limits as the cropped video: MMB_WriteEnFaceVideo takes
        % P(CLimPct) of each movie, so clip the full movies at the cropped
        % limits and ask for P100 (= that limit; the near field saturates)
        CLe = [MMB_SymmetricCLim(Movie, Cfg.WaveCLimPct), MMB_SymmetricCLim(MovieDF, Cfg.WaveCLimPct)];
        Mf  = max(min(Movie_full,   CLe(1)), -CLe(1));
        MDf = max(min(MovieDF_full, CLe(2)), -CLe(2));
        MMB_WriteEnFaceVideo({Mf, MDf}, {'Recovered (near field kept)', 'Recovered + outward DF (near field kept)'}, ...
            xq, yq, Time, [], fullfile(outDir, 'Video_EnFace_noDF_vs_DF_Full.mp4'), 'FrameIdx', 1:nF, ...
            'FrameRate', Cfg.Video.FrameRate, 'Quality', Cfg.Video.Quality, 'CLimPct', 100, ...
            'Units', motionUnits, 'MainTitle', sprintf('SBM %dx%d, %.0f Hz, near field kept |', nA, nB, f_an));
    end
end

%% 14. Save

Info = struct();
Info.Experiment = E;
Info.ScanInfo = ScanInfo;
Info.FileCheck = FileCheck;
Info.SidecarFile = SidecarFile;
Info.SBM = rmfield(SInfo, intersect(fieldnames(SInfo), {'IsAlineBody'}));
Info.f_exc = f_exc;  Info.f_an = f_an;  Info.T_SBM = T_SBM;
Info.MotionUnits = motionUnits;  Info.VelScale = VelScale;
Info.Spectrum = struct('f_equiv_Hz', H.f_equiv_Hz, 'Amp', H.Spectrum, 'SNR_dB', H.SNR_dB, ...
    'PeakBin', H.PeakBin, 'PeakIsExcitation', H.PeakIsExcitation);
Info.Detection = H.Detect;
Info.Band = rmfield(H.Band, intersect(fieldnames(H.Band), {'P', 'P_dB', 'freqs'}));
if isfield(H.Band, 'freqs'), Info.BandSpectrum = struct('freqs', H.Band.freqs, 'P', H.Band.P); end
Info.BulkPhase = BulkPhase;
Info.Source = Src;  Info.Source_mm = Source_mm;  Info.CohMap = CohMap;
Info.nNeedlePx = nnz(needleTip);  Info.nShaftPx = nnz(needleShaft & ~needleTip);
Info.ValidSurfaceFraction = mean(validPix(:));
Info.AcqTime_s = SInfo.AcqWindow_s;
Info.c_ref = c_ref;
Info.Cfg = Cfg;
Info.ProcessingTime_min = toc(tStartAll)/60;
DFInfo.Spec_Before = []; DFInfo.Spec_After = [];      % large, figures already saved

save(fullfile(outDir, 'Results_SBM.mat'), 'Res', 'ResDF', 'Stats', 'StatsDF', 'StatsFar', 'StatsFarDF', ...
    'Ring', 'RingDF', 'Info', 'DFInfo', 'H', ...
    'Border', 'needle', 'SurfMotion', 'FramesRaw', 'U', 'UDF', 'x_mm', 'y_mm', 'Zaxis', 'Time', ...
    'Bmode_IntLog', 'SurfDepth', 'BmodeVol', '-v7.3');
fprintf('    done in %.1f min | medians noDF: %s | DF: %s m/s\n', Info.ProcessingTime_min, ...
    mat2str(round(Stats.Median', 2)), mat2str(round(StatsDF.Median', 2)));
if ~Cfg.ShowFigures, close all; end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function fig = plotSpaceTime(STs, xs, Time, CL, xNeedle, yLine, units, mainTitle, suffix)
% Three stacked space-time panels (measured / recovered / recovered + DF),
% MMB_RedBlue, NaN mid-gray, symmetric limits +-CL, needle x dotted.
fig = figure('Position', [40 40 1500 950], 'Color', 'w');
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
ttl = {'Measured frames (frame order = equivalent time, before compensation)', ...
       ['Recovered (Eq. 59)' suffix], ['Recovered + outward directional filter' suffix]};
for i = 1:3
    ax = nexttile(tl);
    h_ = imagesc(ax, xs{i}, Time, STs{i}.'); set(h_, 'AlphaData', isfinite(STs{i}.'));
    set(ax, 'Color', [0.6 0.6 0.6]); colormap(ax, MMB_RedBlue(255)); clim(ax, [-CL CL]);
    xline(ax, xNeedle, 'k:'); ylabel(ax, 'equivalent time (ms)'); title(ax, ttl{i}); set(ax, 'FontSize', 12);
end
xlabel(ax, sprintf('x (mm) at y = %.2f mm (through the needle)', yLine));
cb = colorbar(ax); cb.Layout.Tile = 'east'; ylabel(cb, units);
title(tl, mainTitle, 'FontWeight', 'bold');
end
