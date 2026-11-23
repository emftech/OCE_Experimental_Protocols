%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% Make_WaveFigures_RedBlue_MMB
%
% Regenerates the wave-propagation displays of every processed experiment
% with the blue-white-red colormap (white = zero motion, MMB_RedBlue):
%   * SpaceTimeMap_noDF_vs_DF (.png/.fig)  - from the saved space-time maps
%   * Video_EnFace_noDF_vs_DF.mp4           - the en face movies are
%     re-interpolated from the saved narrowband traces (Frames, FramesDF)
%     with the same settings used by the 2D estimators (100 x 100 grid,
%     linear Delaunay, Gaussian 0.8 px), near-field disk masked.
% Elastograms and the k-f directional-filter figure keep turbo. No raw data
% is read and no result changes.
% Since 20-10-2026 (later) every display is written twice: cropped (needle
% near field / shaft in gray) and '_Full' (near field kept, SAME color
% limit as the cropped version, so the near field saturates and the far
% field reads exactly as in the cropped version).
%
% Workspace options: OnlyExperiments = {...}; WorkerID = k; NumWorkers = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;
Root = fileparts(mfilename('fullpath'));   % repository root
addpath(genpath(fullfile(Root, 'Codes')));
set(0, 'DefaultFigureVisible', 'off');
if ~exist('WorkerID', 'var'), WorkerID = 1; end
if ~exist('NumWorkers', 'var'), NumWorkers = 1; end
WaveCLimPct = 97;

d = dir(fullfile(Root, 'Results', 'PC_*', 'Results_MMB.mat'));
names = arrayfun(@(x) char(regexp(x.folder, 'PC_[^/\\]+', 'match', 'once')), d, 'UniformOutput', false);
if exist('OnlyExperiments', 'var') && ~isempty(OnlyExperiments)
    names = names(ismember(names, OnlyExperiments));
end
names = names(WorkerID:NumWorkers:end);

for e = 1:numel(names)
    name = names{e};
    tS = tic;
    outDir = fullfile(Root, 'Results', name);
    S = load(fullfile(outDir, 'Results_MMB.mat'), 'Info', 'Frames', 'FramesDF', 'needle', 'SpaceTime', ...
        'SpaceTimeDF', 'Xaxis', 'Time', 'Geom');
    Cfg = S.Info.Cfg; G = S.Geom;
    nPos = G.nPosPerMer; nMer = G.nMer; f_an = S.Info.f_an;

    % Space-time maps. Cropped: needle tip disk + shaft rows NaN (gray),
    % CLim = P97 of the cropped no-DF map. Full: the same data WITHOUT
    % the needle blanking (DF re-run on the full traces), SAME CLim, so the
    % near field saturates and the far field looks identical.
    ttlST = sprintf('Space-time maps - %d meridians, f_{an} = %.0f Hz', nMer, f_an);
    [fig, CLst] = MMB_PlotSpaceTime({S.SpaceTime, S.SpaceTimeDF}, {'Narrowband filtered', 'Narrowband + directional (outward) filter'}, ...
        S.Xaxis, S.Time, S.Info.ScanLength_mm, ttlST, 'CLimPct', WaveCLimPct);
    MMB_SaveFigure(fig, fullfile(outDir, 'SpaceTimeMap_noDF_vs_DF'));
    close(fig);

    FramesFull = S.Frames;                                   % needle rows kept
    FramesDFFull = MMB_DirectionalFilterOutward(FramesFull, G, S.Info.fs, f_an, ...
        'SigmaKappaRel', Cfg.DF.SigmaKappaRel, 'SigmaFreqBins', Cfg.DF.SigmaFreqBins, ...
        'PadFactor', Cfg.DF.PadFactor, 'Source_mm', S.Info.Source_mm);
    STfull   = FramesFull   - mean(FramesFull,   2, 'omitnan');
    STdfFull = FramesDFFull - mean(FramesDFFull, 2, 'omitnan');
    fig = MMB_PlotSpaceTime({STfull, STdfFull}, {'Narrowband filtered (near field kept)', ...
        'Narrowband + directional (outward) filter (near field kept)'}, S.Xaxis, S.Time, S.Info.ScanLength_mm, ...
        [ttlST ' - full, cropped-version limits'], 'CLim', CLst);
    MMB_SaveFigure(fig, fullfile(outDir, 'SpaceTimeMap_noDF_vs_DF_Full'));
    close(fig);

    % En face movies (same interpolation as OCE_PG_ProcessMB)
    FramesIn = S.Frames; FramesIn(S.needle, :) = NaN;
    srcs = {FramesIn, S.FramesDF, FramesFull, FramesDFFull};
    Mv = cell(1, 4);
    for k = 1:4
        [Mv{k}, Interp] = OCE_MB_InterpolateEnFaceFrames(OCE_MB_ReshapeToMeridians(srcs{k}, nPos, nMer), G, ...
            'GridSize', Cfg.PG.GridSize, 'InterpMethod', 'linear', 'Smoothing', true, 'Sigma_px', 0.8, 'Verbose', false);
    end
    [Xq, Yq] = meshgrid(Interp.x_mm, Interp.y_mm);
    Near = hypot(Xq - S.Info.Source_mm(1), Yq - S.Info.Source_mm(2)) < Cfg.NeedleRadius_mm;
    for k = 1:2
        Mv{k}(repmat(Near, 1, 1, size(Mv{k}, 3))) = NaN;
    end
    vArgs = {'FrameIdx', Cfg.Video.Time_ini:size(Mv{1}, 3), 'FrameRate', Cfg.Video.FrameRate, ...
        'Quality', Cfg.Video.Quality, 'MainTitle', sprintf('%d mer, %.0f Hz |', nMer, f_an)};
    CLef = MMB_WriteEnFaceVideo(Mv(1:2), {'No DF', 'With DF (outward)'}, Interp.x_mm, Interp.y_mm, S.Time, G, ...
        fullfile(outDir, 'Video_EnFace_noDF_vs_DF.mp4'), 'CLimPct', WaveCLimPct, vArgs{:});
    MMB_WriteEnFaceVideo(Mv(3:4), {'No DF (near field kept)', 'With DF (near field kept)'}, Interp.x_mm, Interp.y_mm, S.Time, G, ...
        fullfile(outDir, 'Video_EnFace_noDF_vs_DF_Full.mp4'), 'CLim', CLef, vArgs{:});
    fprintf('[%d/%d] %s done in %.1f min\n', e, numel(names), name, toc(tS)/60);
end
