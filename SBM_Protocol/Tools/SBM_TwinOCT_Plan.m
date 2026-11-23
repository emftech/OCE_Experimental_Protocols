function T = SBM_TwinOCT_Plan(FOV_mm, SamplesPerAline)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_Plan
%
% Planning table for SBM on the Twin OCT: for the thesis cases (11/21/41/81
% B-scans) and for slower-galvo variants (N_A + N_turn = 2P, 4P) it lists
% the GUI values, acquisition time, file size and the X-galvo line speed.
% Uses the period / frames / periods in SBM_TwinOCT_Config.m.
%
%   T = SBM_TwinOCT_Plan(12, 2880)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1, FOV_mm = 12; end
if nargin < 2, SamplesPerAline = 2880; end
Cfg0 = SBM_TwinOCT_Config();
Cfg0.ExcMode = "sine";
P = round(Cfg0.T_excitation/Cfg0.Ts_OCT);
NPP = Cfg0.nPeriods*P;
fprintf('\nT_exc = %.0f us, P = %d. Feasible N_F for N_P = %d: %s\n', ...
    Cfg0.T_excitation*1e6, P, Cfg0.nPeriods, mat2str(find(mod(NPP,1:NPP)==0 & NPP./(1:NPP) < P/2)));

geo = [37 20; 94 20; 208 20];          % (N_A, N_turn): S = P, 2P, 4P for P = 57
geo = geo(mod(sum(geo,2), P) == 0, :);
if isempty(geo), geo = [P-20 20; 2*P-20 20]; end
rows = {};
warning('off', 'SBM:BscanPhase');
for g = 1:size(geo,1)
    for nB = [11 21 41 81]
        Cfg = Cfg0;  Cfg.nAscans = geo(g,1);  Cfg.nTurnSamples = geo(g,2);  Cfg.nBscans = nB;
        try
            [~, I] = SBM_TwinOCT_Waveform(Cfg);
        catch ME
            fprintf('  skip N_A=%d N_B=%d: %s\n', Cfg.nAscans, nB, ME.message); continue;
        end
        lineTime = Cfg.nAscans*Cfg.Ts_OCT;
        rows(end+1,:) = {Cfg.nAscans, Cfg.nTurnSamples, nB, Cfg.nAscans*nB, I.R, I.C, ...
            Cfg.nFrames, I.F, round(I.AcqWindow_s*1e3,1), round(FOV_mm/(lineTime*1e3),1), ...
            1/(2*(Cfg.nAscans+Cfg.nTurnSamples)*Cfg.Ts_OCT), ...
            round(2*SamplesPerAline*I.TotalAlines/2^20)}; %#ok<AGROW>
    end
end
warning('on', 'SBM:BscanPhase');
T = cell2table(rows, 'VariableNames', {'N_A','N_turn','N_B','PixelsPerFrame', ...
    'GUI_AlinesPerBframe','GUI_Bframes','GUI_NoCluster','FramePeriod_samples', ...
    'AcqTime_ms','XGalvoSpeed_mm_per_ms','XTriangleFreq_Hz','FileSize_MiB'});
disp(T);
end
