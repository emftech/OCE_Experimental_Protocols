function E = MMB_ParseExperimentName(name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_ParseExperimentName
%
% Decodes the acquisition parameters written in the file name of a
% Multimeridian MB-mode experiment (same columns as the experiment log):
%
%   PC_QH_285us_3c_1Vpp_450_50_100_OCE_8_30_12_12_r1
%   |  |  |     |  |    |   |  |   |   | |  |  |  repetition
%   |  |  |     |  |    |   |  |   |   | |  |  yW [mm]
%   |  |  |     |  |    |   |  |   |   | |  xW [mm]
%   |  |  |     |  |    |   |  |   |   | B-frame signal delay
%   |  |  |     |  |    |   |  |   |   meridians
%   |  |  |     |  |    |   |  |   scan type
%   |  |  |     |  |    |   |  3D number of B-frames (positions/meridian)
%   |  |  |     |  |    |   no. sync samples
%   |  |  |     |  |    A-lines per B-frame (M-mode repetitions)
%   |  |  |     |  Vpp (p = decimal point)
%   |  |  |     number of cycles
%   |  |  excitation period
%   |  temporal pattern (QH)
%   phantom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, name] = fileparts(char(name));
tk = split(string(name), '_');

E = struct();
E.Name        = name;
E.Phantom     = char(tk(1));
E.Pattern     = char(tk(2));
E.PeriodStr   = char(tk(3));
E.Period_us   = str2double(erase(tk(3), 'us'));
E.f_exc       = 1e6/E.Period_us;
E.Cycles      = str2double(erase(tk(4), 'c'));
E.Vpp         = str2double(strrep(erase(tk(5), 'Vpp'), 'p', '.'));
E.Alines      = str2double(tk(6));
E.SyncSamples = str2double(tk(7));
E.Bframes     = str2double(tk(8));
E.ScanType    = char(tk(9));
E.Meridians   = str2double(tk(10));
E.BframeDelay = str2double(tk(11));
E.xW_mm       = str2double(tk(12));
E.yW_mm       = str2double(tk(13));
E.Repetition  = char(tk(14));
E.RepIdx      = str2double(erase(tk(14), 'r'));

end
