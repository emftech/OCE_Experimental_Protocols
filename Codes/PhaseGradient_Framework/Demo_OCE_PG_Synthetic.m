%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% Demo_OCE_PG_Synthetic
%
% Validation of the experimental phase-gradient framework on synthetic
% Multimeridian MB-mode data with known ground truth. The synthetic data
% has EXACTLY the layout produced by the experimental pipeline
% (loaded_phases_Border: [nMer*nPos x Nt], meridians concatenated in
% acquisition order), with the acquisition parameters of the 25-09-2026
% experiments (10 mm meridians, 100 positions, 950 A-lines at 200 kHz,
% 1 kHz excitation).
%
% Cases:
%   A) Cylindrical wave from the pattern center, homogeneous c = 3 m/s
%   B) Non-cylindrical wave: source displaced from the pattern center.
%      The wave vector is not parallel to the meridians, so 1D PG
%      overestimates the speed (c / cos(angle)), while 2D PG does not.
%   C) Cylindrical wave in a concentric two-speed medium:
%      c = 2 m/s for r < 2.5 mm and 4 m/s outside (exact radial phase).
%
% Run from any folder once the framework is on the MATLAB path.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;

if ~exist('SaveDir', 'var')
    SaveDir = '';
end

%% ===================== Acquisition parameters ============================

nPos = 100;             % Bframes_in_3Dscan
nMer = 8;               % No_3Dscans
scanLength_mm = 10;     % Ver_scan_length_mm
fs = 200e3;             % A-line rate [Hz]
Nt = 949;               % time samples after diff() of 950 A-lines
f_exc = 1000;           % [Hz]
omega0 = 2*pi*f_exc;
SNR_dB = 10;            % per-trace SNR
fracNaNTraces = 0.02;   % traces lost (e.g. surface not found)

rng(1);

Geom = OCE_MB_MeridianGeometry(nPos, nMer, scanLength_mm);
t = (0:Nt-1)/fs;

%% ===================== Synthetic cases ===================================

Cases = struct([]);

% A) cylindrical, homogeneous
Cases(1).Name = 'A: cylindrical, c = 3 m/s';
Cases(1).Source_mm = [0 0];
Cases(1).TravelTime = @(r_mm) (r_mm*1e-3)/3;
Cases(1).Speed = @(X,Y) 3*ones(size(X));

% B) non-cylindrical, homogeneous (source displaced)
Cases(2).Name = 'B: displaced source, c = 3 m/s';
Cases(2).Source_mm = [-7 -3];
Cases(2).TravelTime = @(r_mm) (r_mm*1e-3)/3;
Cases(2).Speed = @(X,Y) 3*ones(size(X));

% C) cylindrical, concentric two-speed medium
Rinc = 2.5; c_in = 2; c_out = 4;
Cases(3).Name = 'C: concentric medium, 2 / 4 m/s';
Cases(3).Source_mm = [0 0];
Cases(3).TravelTime = @(r_mm) (min(r_mm,Rinc)*1e-3)/c_in + (max(r_mm-Rinc,0)*1e-3)/c_out;
Cases(3).Speed = @(X,Y) c_in*(hypot(X,Y) < Rinc) + c_out*(hypot(X,Y) >= Rinc);

%% ===================== Run the framework on every case ===================

Summary = table();

for ic = 1:numel(Cases)

    C = Cases(ic);
    fprintf('\n=============== Case %s ===============\n', C.Name);

    % Wavefield at every MB sample, in acquisition order
    r_src = hypot(Geom.X_mm(:) - C.Source_mm(1), Geom.Y_mm(:) - C.Source_mm(2));
    T = C.TravelTime(r_src);
    A = 1 ./ sqrt(r_src + 0.3);
    SpaceTime = A .* sin(omega0*(t - T));

    % Noise and lost traces
    sigma = rms(SpaceTime(:)) / 10^(SNR_dB/20);
    SpaceTime = SpaceTime + sigma*randn(size(SpaceTime));
    lost = rand(size(SpaceTime,1),1) < fracNaNTraces;
    SpaceTime(lost,:) = NaN;

    Results = OCE_PG_ProcessMB(SpaceTime, fs, f_exc, scanLength_mm, nPos, nMer, ...
        'KernelSize_px', 11, 'SpeedLimits', [0.5 20], 'Verbose', true);

    % Ground truth on the en face grid
    [Xq, Yq] = meshgrid(Results.x_mm, Results.y_mm);
    GT = C.Speed(Xq, Yq);

    names = fieldnames(Results.Maps);
    for i = 1:numel(names)
        E = Results.Maps.(names{i});
        v = isfinite(E);
        Ea = 100*abs(E(v) - GT(v))./GT(v);
        row = table(string(C.Name), string(names{i}), median(E(v)), ...
            median(Ea), mean(Ea), 100*nnz(v)/nnz(Results.PG1D.SupportMask), ...
            'VariableNames', {'Case','Method','MedianSpeed','MedianEa_pct','MeanEa_pct','Coverage_pct'});
        Summary = [Summary; row]; %#ok<AGROW>
    end

    if ic == 3
        Inc = hypot(Xq, Yq) < Rinc - 0.6;
        Bkg = hypot(Xq, Yq) > Rinc + 0.6;
        fprintf('\nCase C ROI statistics:\n');
        disp(OCE_PG_ROIStats(Results, Inc, Bkg));
        figs = OCE_PG_PlotResults(Results, 'CLim', [1 6], 'Inc_Mask', hypot(Xq,Yq) < Rinc, ...
            'SaveDir', SaveDir, 'Tag', '_CaseC');
    elseif ~isempty(SaveDir)
        OCE_PG_PlotResults(Results, 'CLim', [1 6], 'SaveDir', SaveDir, ...
            'Tag', sprintf('_Case%c', 'A' + ic - 1));
    end

    CaseResults{ic} = Results; %#ok<SAGROW>
end

fprintf('\n=============== Summary (masked maps) ===============\n');
disp(Summary);
