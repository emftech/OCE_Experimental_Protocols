%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% Make_DepthVideos_RepMeridian_MMB
%
% Depth motion video of ONE representative meridian per experiment
% (Results/<exp>/Video_Depth_RepMeridian.mp4), for the slides: a single
% meridian shows the wave leaving the needle much more clearly than all
% meridians concatenated. The meridian is chosen by
% MMB_SelectRepresentativeMeridian (most reliable pixels, ties -> closest
% to the needle tip). Only that meridian is re-read from the raw file;
% band, masks and B-mode come from Results_MMB.mat (nothing else changes).
%
% Workspace options: OnlyExperiments = {...}; WorkerID = k; NumWorkers = n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;
Root = fileparts(mfilename('fullpath'));   % repository root
addpath(genpath(fullfile(Root, 'Codes')));
set(0, 'DefaultFigureVisible', 'off');
if ~exist('WorkerID', 'var'), WorkerID = 1; end
if ~exist('NumWorkers', 'var'), NumWorkers = 1; end

d = dir(fullfile(Root, 'Results', 'PC_*', 'Results_MMB.mat'));
names = arrayfun(@(x) char(regexp(x.folder, 'PC_[^/\\]+', 'match', 'once')), d, 'UniformOutput', false);
if exist('OnlyExperiments', 'var') && ~isempty(OnlyExperiments)
    names = names(ismember(names, OnlyExperiments));
end
names = names(WorkerID:NumWorkers:end);
RepTable = table();

for e = 1:numel(names)
    name = names{e};
    tS = tic;
    outDir = fullfile(Root, 'Results', name);
    S = load(fullfile(outDir, 'Results_MMB.mat'), 'Info', 'Bmode_IntLog', 'VideoMask', 'Zaxis', 'Time', 'Geom');
    Cfg = S.Info.Cfg; SI = S.Info.ScanInfo;
    nS = SI.samples_in_Aline; nA = SI.Alines_in_Bframe; nPos = SI.Bframes_in_3Dscan; nMer = SI.No_3Dscans;
    nTot = nPos*nMer;

    [b, RI] = MMB_SelectRepresentativeMeridian(S.VideoMask, S.Bmode_IntLog, S.Geom, S.Info.Source_mm, 90);
    posIdx = (b-1)*nPos + (1:nPos);

    hit = dir(fullfile(Root, 'Data', '**', [name '.bin']));
    binFile = fullfile(hit(1).folder, hit(1).name);
    [~, hb] = MMB_ReadBinHeader(binFile);
    m = memmapfile(binFile, 'Offset', hb, 'Format', {'uint16', [nS nA nTot], 'raw'});
    tIdx = Cfg.Cut_Time(1):min(Cfg.Cut_Time(2), nA);
    zIdx = Cfg.Cut_Depth(1):min(Cfg.Cut_Depth(2), nS/2);
    nT = numel(tIdx); hannK = hann(nS);
    dc = zeros(nS, nTot);
    for p = 1:nTot
        dc(:,p) = median(double(m.Data.raw(:,:,p)), 2);
    end
    dc = median(dc, 2);

    C = MMB_ComputeCplx(m, posIdx, dc, hannK, tIdx, zIdx);
    v = MMB_GaussianBandpass(single(diff(unwrap(-angle(C), [], 3), 1, 3))*S.Info.VelScale_mmps_per_rad, ...
        S.Info.fs, S.Info.f_an, S.Info.Band.BW_Hz, 3);
    VelMer = permute(v - mean(v, 3), [2 1 3]);              % [z x pos x t]
    clear C v

    MMB_WriteRepMeridianDepthVideo(VelMer, S.Bmode_IntLog(:, posIdx), S.VideoMask(:, posIdx), S.Geom, b, ...
        S.Info.Source_mm, S.Zaxis, S.Time, fullfile(outDir, 'Video_Depth_RepMeridian.mp4'), ...
        'FrameIdx', Cfg.Video.Time_ini:nT-1, 'FrameRate', Cfg.Video.FrameRate, 'Quality', 60, ...
        'CLimPct', 90, 'MinBmode_dB', 90, 'IThresh', Cfg.i_thresh, ...
        'MainTitle', sprintf('%d mer, meridian %d, %.0f Hz |', nMer, b, S.Info.f_an));

    RepMeridian = RI; save(fullfile(outDir, 'RepMeridian.mat'), 'RepMeridian');
    RepTable = [RepTable; table(string(name), b, RI.Angle_deg(b), RI.DistToTip_mm(b), RI.Score(b)/max(RI.Score), ...
        'VariableNames', {'Experiment', 'Meridian', 'Angle_deg', 'DistToTip_mm', 'RelScore'})]; %#ok<AGROW>
    fprintf('[%d/%d] %s: meridian %d (%.1f deg, %.2f mm from tip) in %.1f min\n', e, numel(names), name, ...
        b, RI.Angle_deg(b), RI.DistToTip_mm(b), toc(tS)/60);
end
disp(RepTable);
