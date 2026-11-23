%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% Make_DepthVideos_CroppedAndFull_MMB
%
% Depth motion videos of every processed experiment, in two versions:
%   * cropped (as before): overlay only on tissue below the surface with
%     the needle columns removed; CL = P90 of |v| over reliable pixels
%       Video_Depth_AllMeridians_Masked.mp4, Video_Depth_RepMeridian.mp4
%   * '_Full': the same data with the needle / near-field columns KEPT
%     (overlay on all tissue below the surface) and the SAME CL as the
%     cropped version, so the near field saturates and the far field looks
%     exactly as in the cropped video
%       Video_Depth_AllMeridians_Masked_Full.mp4, Video_Depth_RepMeridian_Full.mp4
% The raw .bin is read once per experiment (depth velocity); band, masks,
% border and B-mode come from Results_MMB.mat, the representative meridian
% from RepMeridian.mat (same choice as before).
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
    S = load(fullfile(outDir, 'Results_MMB.mat'), 'Info', 'Bmode_IntLog', 'Bmode_Mask', 'VideoMask', 'Border', ...
        'Xaxis', 'Zaxis', 'Time', 'Geom');
    Cfg = S.Info.Cfg; SI = S.Info.ScanInfo;
    nS = SI.samples_in_Aline; nA = SI.Alines_in_Bframe; nPos = SI.Bframes_in_3Dscan; nMer = SI.No_3Dscans;
    nTot = nPos*nMer;

    % Depth velocity, all meridians (same steps as section 10 of the main script)
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
    Vel = zeros(nZ, nTot, nT-1, 'single');
    for b = 1:nMer
        posIdx = (b-1)*nPos + (1:nPos);
        C = MMB_ComputeCplx(m, posIdx, dc, hannK, tIdx, zIdx);
        v = MMB_GaussianBandpass(single(diff(unwrap(-angle(C), [], 3), 1, 3))*S.Info.VelScale_mmps_per_rad, ...
            S.Info.fs, S.Info.f_an, S.Info.Band.BW_Hz, 3);
        Vel(:, posIdx, :) = permute(v - mean(v, 3), [2 1 3]);
    end
    clear C v

    % Full mask: tissue below the detected surface, needle columns KEPT
    MaskFull = S.Bmode_Mask;
    for c = 1:nTot
        bI = S.Border.Idx(c);
        if isfinite(bI)
            MaskFull(1:max(round(bI)-1, 0), c) = false;
        else
            MaskFull(:, c) = false;
        end
    end

    common = {'FrameIdx', Cfg.Video.Time_ini:nT-1, 'FrameRate', Cfg.Video.FrameRate, 'MinBmode_dB', 90, ...
        'IThresh', Cfg.i_thresh};

    % All meridians: cropped, then full with the cropped CL
    ttl = sprintf('%d mer, %.0f Hz |', nMer, S.Info.f_an);
    CLall = MMB_WriteDepthVideo(Vel, S.Bmode_IntLog, S.VideoMask, S.Xaxis, S.Zaxis, S.Time, ...
        fullfile(outDir, 'Video_Depth_AllMeridians_Masked.mp4'), common{:}, 'Quality', 50, 'CLimPct', 90, ...
        'MeridianLength_mm', S.Info.ScanLength_mm, 'MainTitle', ttl);
    MMB_WriteDepthVideo(Vel, S.Bmode_IntLog, MaskFull, S.Xaxis, S.Zaxis, S.Time, ...
        fullfile(outDir, 'Video_Depth_AllMeridians_Masked_Full.mp4'), common{:}, 'Quality', 50, 'CLim', CLall, ...
        'MeridianLength_mm', S.Info.ScanLength_mm, 'MainTitle', [ttl ' near field kept |']);

    % Representative meridian: cropped, then full with the cropped CL
    RM = load(fullfile(outDir, 'RepMeridian.mat'));
    b = RM.RepMeridian.Meridian;
    posIdx = (b-1)*nPos + (1:nPos);
    ttl = sprintf('%d mer, meridian %d, %.0f Hz |', nMer, b, S.Info.f_an);
    CLrep = MMB_WriteRepMeridianDepthVideo(Vel(:, posIdx, :), S.Bmode_IntLog(:, posIdx), S.VideoMask(:, posIdx), ...
        S.Geom, b, S.Info.Source_mm, S.Zaxis, S.Time, fullfile(outDir, 'Video_Depth_RepMeridian.mp4'), ...
        common{:}, 'Quality', 60, 'CLimPct', 90, 'MainTitle', ttl);
    MMB_WriteRepMeridianDepthVideo(Vel(:, posIdx, :), S.Bmode_IntLog(:, posIdx), MaskFull(:, posIdx), ...
        S.Geom, b, S.Info.Source_mm, S.Zaxis, S.Time, fullfile(outDir, 'Video_Depth_RepMeridian_Full.mp4'), ...
        common{:}, 'Quality', 60, 'CLim', CLrep, 'MainTitle', [ttl ' near field kept |']);

    clear Vel
    fprintf('[%d/%d] %s done in %.1f min (CL all %.3f, rep %.3f mm/s)\n', e, numel(names), name, toc(tS)/60, CLall, CLrep);
end
