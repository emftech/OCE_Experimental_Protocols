function [waveform, TotalClusterScans, Info] = SBM_TwinOCT_LabVIEWNode( ...
    horizontalMin, horizontalMax, verticalMin, verticalMax, ...
    UI_AlinesPerBframe, UI_Bframes, UI_NoCluster, UI_SyncSamples, UI_BframeDelay, ...
    sidecarDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_LabVIEWNode
%
% Entry point called from the MATLAB script node of
% "Waveform SubVIs/Galvo 3D cluster driver_SBM.vi" (see
% ScriptNode_GalvoDriver_SBM.txt). It
%   1. reads the protocol from SBM_TwinOCT_Config.m,
%   2. takes the galvo voltages and the B-frame Signal Delay from the GUI,
%   3. checks that the Alazar settings typed in the GUI match the protocol
%      (otherwise the acquisition would record the wrong A-lines) and stops
%      with a message listing the values to type,
%   4. returns the AO waveform (rows X, Y, B-trigger[, Excitation]) and the
%      total number of A-lines (TotCluster_Scans),
%   5. saves a sidecar .mat (Cfg + Info + waveform) for the processing.
%
% GUI value          script-node terminal   SBM meaning
% A-lines per B-frame  M_repetitions        R  = A-lines per B-trigger
% B-frames / Vert.     nAscans              C  = B-triggers per SBM frame
% No. Cluster          nBscans              N_F = SBM frames
% No. Sync Samples     waitTimeAscans       N_turn
% B-frame Signal Delay Bframesignaldelay    galvo-lag delay of the trigger
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cfg = SBM_TwinOCT_Config();
Cfg.horizontalMin_V = double(horizontalMin);
Cfg.horizontalMax_V = double(horizontalMax);
Cfg.verticalMin_V   = double(verticalMin);
Cfg.verticalMax_V   = double(verticalMax);
Cfg.TrigDelay       = round(double(UI_BframeDelay));

[waveform, Info] = SBM_TwinOCT_Waveform(Cfg);
L = Info.LabVIEW;

got  = round(double([UI_AlinesPerBframe UI_Bframes UI_NoCluster UI_SyncSamples]));
want = [L.A_lines_per_B_frame L.B_frames L.No_Cluster L.No_Sync_Samples];
if ~isequal(got, want)
    error('SBM:GUI', ['SBM GUI mismatch. Set in Scan Parameters:\n' ...
        '  A-lines per B-frame = %d (now %d)\n  B-frames = %d (now %d)\n' ...
        '  No. Cluster = %d (now %d)\n  No. Sync Samples = %d (now %d)\n' ...
        'or edit SBM_TwinOCT_Config.m.'], ...
        want(1), got(1), want(2), got(2), want(3), got(3), want(4), got(4));
end
TotalClusterScans = Info.TotalAlines;

if nargin >= 10 && strlength(string(sidecarDir)) > 0
    sidecarDir = char(sidecarDir);
    if ~exist(sidecarDir, 'dir'), mkdir(sidecarDir); end
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    save(fullfile(sidecarDir, ['SBM_Waveform_' stamp '.mat']), 'Cfg', 'Info', 'waveform', '-v7.3');
    save(fullfile(sidecarDir, 'SBM_Waveform_LAST.mat'), 'Cfg', 'Info', 'waveform', '-v7.3');
end
end
