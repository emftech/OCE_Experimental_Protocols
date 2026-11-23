%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% Check_ExperimentalChain_On_Simulations_MB
%
% Apples-to-apples check: runs the EXPERIMENTAL multimeridian MB-mode chain
% of Main_Processing_SingleExp_MultimerMB_EMF.m (sections 6-12) on the
% simulated traces of the "Acquisition protocol simulations" session
% (Acquisition Protocols Simulator/Results/Experiment_Matched_18x18mm/...).
%
%   narrowband Gaussian filter at f_exc (-3 dB width, <= 600 Hz)
%   -> near-field disk around the (known) source blanked
%   -> outward directional filter (MMB_DirectionalFilterOutward)
%   -> OCE_PG_ProcessMB (11 px, 100 x 100, R2 >= 0.7, mu = 0.1 median|k|, exact DTFT,
%      neumann TV), without and with DF, plus a no-DF run with the thesis
%      phasor convention (PhasorMethod fftbin, nearest bin of a 1024-pt FFT)
%   -> MMB_MapStats (Ea vs c_ref = 10.54 m/s, Ep) + radial speed profile
%
% Blanking radius: max(1.2 mm, 1.2*SourceFWHM_mm) by default (4.25 mm for
% the 3.54 mm needle, as suggested by the simulation session).
%
% Output: <SimRoot>/MB_ExperimentalChain_Check/ (table + radial profiles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;
Root = fileparts(mfilename('fullpath'));   % repository root
SimRoot = fullfile(Root, '..', 'SBM_Protocol_OCE', 'Results', 'Experiment_Matched_18x18mm');   % simulator repo cloned next to this one
addpath(genpath(fullfile(Root, 'Codes')));
set(0, 'DefaultFigureVisible', 'off');

c_ref = 0.9553*sqrt(0.365e6/(3*1000));
if ~exist('BlankRadiusOverride_mm', 'var'), BlankRadiusOverride_mm = []; end
OutDir = fullfile(SimRoot, 'MB_ExperimentalChain_Check');
if ~exist(OutDir, 'dir'), mkdir(OutDir); end

files = dir(fullfile(SimRoot, 'FWHM_*', 'Multimeridian_MB', '**', 'MB*Mer_Traces_*us.mat'));
fprintf('%d simulated trace files found.\n', numel(files));
methods = {'1D PG', '1D PG-TV', '2D PG', '2D PG-TV'};
names = {'PG1D', 'PG1D_TV', 'PG2D', 'PG2D_TV'};
T = table(); Prof = struct();

for k = 1:numel(files)
    F = fullfile(files(k).folder, files(k).name);
    L = load(F, 'Traces', 'Info');
    I = L.Info; V = double(L.Traces);
    fs = I.fs; f_exc = I.f_exc; Lmm = I.L_mm; nPos = I.N_pos; nMer = I.N_mer;
    src = [0 0]; if isfield(I, 'Source_mm'), src = I.Source_mm; end
    fwhm = NaN; if isfield(I, 'SourceFWHM_mm'), fwhm = I.SourceFWHM_mm; end
    R = BlankRadiusOverride_mm;
    if isempty(R), R = max(1.2, 1.2*fwhm); if ~isfinite(R), R = 1.2; end, end
    caseName = char(regexp(F, 'FWHM_[^/]+', 'match', 'once'));

    % Narrowband selection and zero-phase Gaussian filter (pass 2: without the near field)
    G = OCE_MB_MeridianGeometry(nPos, nMer, Lmm, 'ReverseAlternate', true);
    near = hypot(G.X_mm(:) - src(1), G.Y_mm(:) - src(2)) < R;
    Band = MMB_SelectNarrowBand(V(~near, :), fs, f_exc, 'SearchRel', 0.2, 'MaxBW_Hz', 600, 'MinProminence_dB', 3);
    Fr = MMB_GaussianBandpass(V, fs, f_exc, Band.BW_Hz, 2);
    Fr(near, :) = NaN;

    % Directional filter and estimators
    [FrDF, DFI] = MMB_DirectionalFilterOutward(Fr, G, fs, f_exc, 'SigmaKappaRel', 0.5, 'SigmaFreqBins', 2, ...
        'PadFactor', 4, 'Source_mm', src);
    args = {'KernelSize_px', 11, 'GridSize', 100, 'R2Threshold', 0.7, 'SpeedLimits', [0.5 40], ...
        'ExcludeCenterRadius_mm', 0, 'MuTVRel', 0.1, 'TVIterations', 500, 'Verbose', false};
    warning('off', 'all');
    Res = {OCE_PG_ProcessMB(Fr, fs, f_exc, Lmm, nPos, nMer, args{:}), ...
           OCE_PG_ProcessMB(FrDF, fs, f_exc, Lmm, nPos, nMer, args{:}), ...
           OCE_PG_ProcessMB(Fr, fs, f_exc, Lmm, nPos, nMer, args{:}, 'PhasorMethod', 'fftbin')};
    warning('on', 'all');
    [Xq, Yq] = meshgrid(Res{1}.x_mm, Res{1}.y_mm);
    rq = hypot(Xq - src(1), Yq - src(2));

    % noDF / DF: exact DTFT at f_exc.  noDF_fftbin: phase at the nearest bin
    % of a 1024-point FFT (thesis convention), speed still 2*pi*f_exc/|k|
    tags = {'noDF', 'DF', 'noDF_fftbin'};
    for d = 1:3
        for i = 1:4
            Res{d}.Maps.(names{i})(rq < R) = NaN;
        end
        St = MMB_MapStats(Res{d}.Maps, c_ref);
        for i = 1:4
            T = [T; table(string(caseName), nMer, 1e6/f_exc, f_exc, string(tags{d}), string(methods{i}), R, ...
                Band.BW3dB_Hz, St.Median(i), St.Mean(i), St.Ep(i), St.Ea(i), St.ValidPx(i), median(DFI.OutwardEnergyFraction), ...
                'VariableNames', {'Case', 'Meridians', 'Period_us', 'f_exc_Hz', 'Filter', 'Method', 'Blank_mm', ...
                'BW3dB_Hz', 'Median_mps', 'Mean_mps', 'Ep_pct', 'Ea_pct', 'ValidPx', 'OutwardFraction'})]; %#ok<SAGROW>
        end
        % Radial profile of 2D PG-TV (median per 0.5 mm ring)
        E = Res{d}.Maps.PG2D_TV; rr = 0:0.5:Lmm/2;
        pr = arrayfun(@(a) median(E(rq >= a & rq < a + 0.5 & isfinite(E))), rr(1:end-1));
        key = matlab.lang.makeValidName(sprintf('%s_%dmer_%dus_%s', caseName, nMer, round(1e6/f_exc), tags{d}));
        Prof.(key) = struct('r_mm', rr(1:end-1) + 0.25, 'c_mps', pr);
    end
    fprintf('[%d/%d] %s %d mer %.0f us: 2D PG-TV noDF %.2f / DF %.2f / fftbin %.2f m/s (blank %.2f mm)\n', k, numel(files), caseName, ...
        nMer, 1e6/f_exc, T.Median_mps(end-8), T.Median_mps(end-4), T.Median_mps(end), R);
end

writetable(T, fullfile(OutDir, 'ExperimentalChain_on_Simulations_MB.csv'));
writetable(T, fullfile(OutDir, 'ExperimentalChain_on_Simulations_MB.xlsx'));
save(fullfile(OutDir, 'ExperimentalChain_on_Simulations_MB.mat'), 'T', 'Prof', 'c_ref');
fprintf('Saved to %s\n', OutDir);
