%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% Make_DepthVideos_RedBlue_MMB
%
% Regenerates ONLY the masked all-meridian depth videos of the 06-10-2026
% experiments with the new overlay of MMB_WriteDepthVideo (blue-white-red,
% transparent at v = 0, limits = P90 of |v| over pixels with B-mode >= 90
% dB). The elastograms and every other result are unchanged, so the rest
% of the pipeline is not re-run: the band, surface, needle mask and B-mode
% are read from Results/<exp>/Results_MMB.mat and only the depth-resolved
% velocity is recomputed from the raw .bin (same steps as section 10 of
% Main_Processing_SingleExp_MultimerMB_EMF.m).
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
mer = cellfun(@(n) str2double(regexp(n, '(?<=_OCE_)\d+', 'match', 'once')), names);
[~, o] = sort(mer, 'descend');
names = names(o);
names = names(WorkerID:NumWorkers:end);

for e = 1:numel(names)
    name = names{e};
    tS = tic;
    outDir = fullfile(Root, 'Results', name);
    S = load(fullfile(outDir, 'Results_MMB.mat'), 'Info', 'Bmode_IntLog', 'VideoMask', 'Xaxis', 'Zaxis', 'Time');
    Cfg = S.Info.Cfg; SI = S.Info.ScanInfo;
    nS = SI.samples_in_Aline; nA = SI.Alines_in_Bframe; nPos = SI.Bframes_in_3Dscan; nMer = SI.No_3Dscans;
    nTot = nPos*nMer;
    hit = dir(fullfile(Root, 'Data', '**', [name '.bin']));
    binFile = fullfile(hit(1).folder, hit(1).name);
    [~, hb] = MMB_ReadBinHeader(binFile);
    m = memmapfile(binFile, 'Offset', hb, 'Format', {'uint16', [nS nA nTot], 'raw'});
    tIdx = Cfg.Cut_Time(1):min(Cfg.Cut_Time(2), nA);
    zIdx = Cfg.Cut_Depth(1):min(Cfg.Cut_Depth(2), nS/2);
    nT = numel(tIdx); nZ = numel(zIdx); hannK = hann(nS);

    dc = zeros(nS, nTot);
    for p = 1:nTot
        dc(:,p) = median(double(m.Data.raw(:,:,p)), 2);
    end
    dc = median(dc, 2);

    fs = S.Info.fs; f_an = S.Info.f_an; BW = S.Info.Band.BW_Hz; VS = S.Info.VelScale_mmps_per_rad;
    Vel = zeros(nZ, nTot, nT-1, 'single');
    for b = 1:nMer
        posIdx = (b-1)*nPos + (1:nPos);
        C = MMB_ComputeCplx(m, posIdx, dc, hannK, tIdx, zIdx);
        v = MMB_GaussianBandpass(single(diff(unwrap(-angle(C), [], 3), 1, 3))*VS, fs, f_an, BW, 3);
        Vel(:, posIdx, :) = permute(v - mean(v, 3), [2 1 3]);
    end
    clear C v

    MMB_WriteDepthVideo(Vel, S.Bmode_IntLog, S.VideoMask, S.Xaxis, S.Zaxis, S.Time, ...
        fullfile(outDir, 'Video_Depth_AllMeridians_Masked.mp4'), ...
        'FrameIdx', Cfg.Video.Time_ini:nT-1, 'FrameRate', Cfg.Video.FrameRate, 'Quality', 50, ...
        'CLimPct', 90, 'MinBmode_dB', 90, 'IThresh', Cfg.i_thresh, ...
        'MeridianLength_mm', S.Info.ScanLength_mm, 'MainTitle', sprintf('%d mer, %.0f Hz |', nMer, f_an));
    clear Vel
    fprintf('[%d/%d] %s done in %.1f min\n', e, numel(names), name, toc(tS)/60);
end
