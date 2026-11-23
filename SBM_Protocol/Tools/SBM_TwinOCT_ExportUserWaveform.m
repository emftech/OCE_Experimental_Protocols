function [file, Info] = SBM_TwinOCT_ExportUserWaveform(Cfg, outFile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_ExportUserWaveform
%
% Fallback route that does not need the MATLAB script node: writes the SBM
% AO waveform as a tab-delimited text file (one row per AO channel, one
% column per A-line clock tick, "%.4f") for the master VI's "User Defined
% Waveform" / "Waveform Path" option (Read Delimited Spreadsheet.vi in
% LabVIEW_Setup_Mirrors_and_Btrig.vi), plus a sidecar .mat with Cfg/Info.
%
% Cfg must hold the galvo VOLTS of your scan (Horizontal/Vertical Min/Max
% (V) shown in the GUI for the chosen scan length).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[waveform, Info] = SBM_TwinOCT_Waveform(Cfg);
[p, n] = fileparts(outFile);
if isempty(p), p = pwd; end
file = fullfile(p, [n '.txt']);
writematrix(round(waveform, 4), file, 'Delimiter', 'tab', 'FileType', 'text');
save(fullfile(p, [n '_Sidecar.mat']), 'Cfg', 'Info', 'waveform', '-v7.3');
fprintf('Wrote %s (%d channels x %d samples)\n', file, size(waveform,1), size(waveform,2));
disp(Info.LabVIEW);
end
