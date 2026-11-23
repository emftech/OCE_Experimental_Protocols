function E = SBM_ParseExperimentName(name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_ParseExperimentName
%
% Decodes an SBM file name built with the convention of
% Experiment_Log_SBM_*.xlsx (same CONCATENATE order as the MB log):
%
%   Phantom_Type_Period_Cycles_Vpp_NA_NT_NB_SBM_NF_Delay_xW_yW_rep
%   PC_H_285us_1c_1Vpp_37_20_21_SBM_57_30_12_12_r1
%
%   NA = A-lines per B-scan, NT = turn samples (No. Sync Samples),
%   NB = B-scans per frame, NF = SBM frames (No. Cluster).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, name] = fileparts(char(name));
t = strsplit(name, '_');
if numel(t) < 14 || ~strcmpi(t{9}, 'SBM')
    error('SBM:Name', ['"%s" does not follow the SBM convention ' ...
        'Phantom_Type_Period_Cycles_Vpp_NA_NT_NB_SBM_NF_Delay_xW_yW_rep.'], name);
end
E = struct();
E.Name       = name;
E.Phantom    = t{1};
E.Pattern    = t{2};
E.Period_us  = str2double(erase(t{3}, 'us'));
E.f_exc      = 1/(E.Period_us*1e-6);
E.Cycles     = t{4};
E.Vpp        = str2double(erase(t{5}, 'Vpp'));
E.NA         = str2double(t{6});
E.NT         = str2double(t{7});
E.NB         = str2double(t{8});
E.ScanType   = t{9};
E.NF         = str2double(t{10});
E.Delay      = str2double(t{11});
E.xW_mm      = str2double(t{12});
E.yW_mm      = str2double(t{13});
E.Repetition = t{14};
E.RepIdx     = str2double(erase(t{14}, 'r'));
end
