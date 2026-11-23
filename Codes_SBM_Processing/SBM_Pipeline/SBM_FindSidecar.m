function [Info, src] = SBM_FindSidecar(binFile, E, ScanInfo, sidecarDirs, protocolDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_FindSidecar
%
% Finds the acquisition sidecar (Info of SBM_TwinOCT_Waveform, saved by the
% LabVIEW script node) that belongs to an SBM .bin file. Search order:
%   1. <bin folder>/<name>_Sidecar.mat or <name>.mat
%   2. SBM_Waveform_yyyyMMdd_HHmmss.mat in sidecarDirs whose protocol matches
%      the file (R, C, N_F from the header; N_A, N_B, N_turn, T_exc from the
%      name), not later than the acquisition start written in the header
%      + 1 min, the CLOSEST in time first (the LabVIEW node stamps the
%      sidecar within ~1 s of the header time)
%   3. SBM_Waveform_LAST.mat in sidecarDirs if it matches
%   4. rebuilt with SBM_TwinOCT_Waveform (protocolDir on the path) from the
%      name: warns, valid only if the config defaults were unchanged
% The chosen Info is always checked against the header and the name.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4, sidecarDirs = {}; end
if nargin < 5, protocolDir = ''; end
if ischar(sidecarDirs) || isstring(sidecarDirs), sidecarDirs = cellstr(sidecarDirs); end
[fdir, fname] = fileparts(binFile);

% 1. next to the .bin
for c = {fullfile(fdir, [fname '_Sidecar.mat']), fullfile(fdir, [fname '.mat'])}
    if isfile(c{1})
        Info = loadInfo(c{1});
        if matches(Info), src = c{1}; return; end
        warning('SBM:Sidecar', '%s does not match %s; ignored.', c{1}, fname);
    end
end

% 2./3. sidecar folders
cand = struct('file', {}, 't', {});
for d = [{fdir}, sidecarDirs(:)']
    if isempty(d{1}) || ~isfolder(d{1}), continue; end
    L = dir(fullfile(d{1}, 'SBM_Waveform_*.mat'));
    for i = 1:numel(L)
        tk = regexp(L(i).name, 'SBM_Waveform_(\d{8}_\d{6})\.mat', 'tokens', 'once');
        if isempty(tk), t = NaT; else, t = datetime(tk{1}, 'InputFormat', 'yyyyMMdd_HHmmss'); end
        cand(end+1) = struct('file', fullfile(L(i).folder, L(i).name), 't', t); %#ok<AGROW>
    end
end
if ~isempty(cand)
    stamped = cand(~isnat([cand.t]));
    if ~isempty(stamped)
        if ~isnat(ScanInfo.Datetime)
            stamped = stamped([stamped.t] <= ScanInfo.Datetime + minutes(1));
            [~, o] = sort(abs([stamped.t] - ScanInfo.Datetime), 'ascend');
        else
            [~, o] = sort([stamped.t], 'descend');
        end
        for i = o
            Info = loadInfo(stamped(i).file);
            if matches(Info), src = stamped(i).file; return; end
        end
    end
    last = cand(isnat([cand.t]) & endsWith({cand.file}, 'SBM_Waveform_LAST.mat'));
    for i = 1:numel(last)
        Info = loadInfo(last(i).file);
        if matches(Info), src = last(i).file; return; end
    end
end

% 4. rebuild
if ~isempty(protocolDir), addpath(protocolDir); end
if exist('SBM_TwinOCT_Waveform', 'file') ~= 2
    error('SBM:NoSidecar', ['No matching sidecar for %s and SBM_TwinOCT_Waveform is not ' ...
        'on the path to rebuild it (set Cfg.ProtocolDir).'], fname);
end
warning('SBM:Rebuilt', ['No matching sidecar for %s: Info rebuilt from the file name and ' ...
    'the SBM_TwinOCT_Config defaults (valid only if those were not changed).'], fname);
Cfg = SBM_TwinOCT_Config();
Cfg.T_excitation = E.Period_us*1e-6;  Cfg.nAscans = E.NA;  Cfg.nTurnSamples = E.NT;
Cfg.nBscans = E.NB;  Cfg.nFrames = E.NF;  Cfg.TrigDelay = E.Delay;
Cfg.horizontalMin_V = -E.xW_mm/2;  Cfg.horizontalMax_V = E.xW_mm/2;
Cfg.verticalMin_V   = -E.yW_mm/2;  Cfg.verticalMax_V   =  E.yW_mm/2;
[~, Info] = SBM_TwinOCT_Waveform(Cfg);
if ~matches(Info)
    error('SBM:Rebuilt', 'Rebuilt Info does not match the header of %s.', fname);
end
src = 'rebuilt from the file name';

    function I = loadInfo(f)
        S = load(f, 'Info');
        if ~isfield(S, 'Info'), I = struct('R', NaN); else, I = S.Info; end
    end
    function ok = matches(I)
        ok = isfield(I, 'Cfg') && I.R == ScanInfo.Alines_in_Bframe && ...
             I.C == ScanInfo.Bframes_in_3Dscan && I.Cfg.nFrames == ScanInfo.No_3Dscans && ...
             I.Cfg.nAscans == E.NA && I.Cfg.nBscans == E.NB && ...
             I.Cfg.nTurnSamples == E.NT && abs(I.Cfg.T_excitation*1e6 - E.Period_us) < 1e-6;
    end
end
