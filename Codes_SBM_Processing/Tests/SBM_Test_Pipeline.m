%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_Test_Pipeline
%
% End-to-end test of Codes_SBM_Processing on synthetic data with a known
% ground truth (c = 10.54 m/s, needle at (0.70, 1.10) mm, 30 % boundary
% reflection): writes 7 SBM .bin files in the LabVIEW format into tempdir
% (tests r1-r3, B-scan sweep 11/41, frequency sweep 220/300 us; sidecars
% next to the .bin, time-stamped and LAST) and runs the batch script on
% them. Expected (19-10-2026 run): medians 10.55-10.63 m/s for 21+ B-scans,
% needle recovered at (0.70, 1.10) mm, all sidecars matched, summary in
% <tempdir>/SBM_Experiments_Synthetic/Results/Summary_SBM.
% Needs the SBM_Protocol folder (SBM_TwinOCT_Waveform) to build the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root  = fileparts(fileparts(mfilename('fullpath')));            % Codes_SBM_Processing
addpath(genpath(root));
proto = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'SBM_Protocol');
ExpDirT = fullfile(tempdir, 'SBM_Experiments_Synthetic');
if exist(ExpDirT, 'dir'), rmdir(ExpDirT, 's'); end
Exps = struct('Period_us', {285, 285, 285, 285, 285, 220, 300}, 'NB', {21, 21, 21, 11, 41, 21, 21}, ...
              'Rep', {1, 2, 3, 1, 1, 1, 1}, ...
              'Sidecar', {'next', 'stamped', 'stamped', 'last', 'next', 'stamped', 'next'});
SBM_MakeSyntheticDataset(ExpDirT, Exps, proto);
CfgOverride = struct('ExpDir', ExpDirT, 'DataDir', fullfile(ExpDirT, 'Data'), ...
    'SidecarDirs', {{fullfile(ExpDirT, 'Sidecars'), fullfile(ExpDirT, 'Data')}}, ...
    'ResultsDir', fullfile(ExpDirT, 'Results'), ...
    'DepthVideo', struct('CLimPct', 90, 'MinBmode_dB', 40), ...  % synthetic B-mode is ~50-57 dB
    'NeedleRadius_mm', 1.2);                                     % point-like synthetic source (thin needle)
Main_Processing_Batch_SBM_EMF;
