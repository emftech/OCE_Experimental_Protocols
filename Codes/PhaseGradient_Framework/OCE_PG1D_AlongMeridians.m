function R = OCE_PG1D_AlongMeridians(Meridian_Frames, Geom, fs, f_exc, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_PG1D_AlongMeridians
%
% 1D Phase-Gradient (1D PG) and TV-regularized 1D Phase-Gradient (1D PG-TV)
% speed estimation for EXPERIMENTAL Multimeridian MB-mode OCE data.
%
% Experimental counterpart of PhaseGradient_1D_AlongMeridians_MBMode_v3.m
% and PhaseGradient_1D_AlongMeridians_Regularized_TV.m (simulation study):
%
%   1) Harmonic phasor of every M-mode trace at f_exc (OCE_HarmonicPhasor).
%   2) Local 1D phase-line fit along each meridian (sliding window of L px):
%
%          theta(s) ~ k s + theta0,    c = omega0 / |k|
%
%   3) Projection of the sparse meridian estimates to the en face plane.
%   4) 1D PG  : interpolation of the sparse speed map  c(x,y).
%      1D PG-TV: interpolation of the sparse |k| map, then 2D TV on K:
%
%          K_TV = argmin_K 0.5 ||K - K_raw||^2 + mu TV(K),  c_TV = w0 / K_TV
%
%   5) Quality masking with the local fit R^2 (thesis: R^2 > 0.7).
%
% Maps are returned on a regular en face grid in MESHGRID convention:
% rows = y, columns = x. Display with imagesc(R.x_mm, R.y_mm, map); axis xy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   Meridian_Frames : [nPosPerMer x Nt x nMer] motion traces (e.g. particle
%                     velocity from MotionAnalysisBorder, reshaped with
%                     OCE_MB_ReshapeToMeridians)
%   Geom            : struct from OCE_MB_MeridianGeometry
%   fs              : A-line rate [Hz]
%   f_exc           : excitation frequency [Hz]
%
% Optional name-value parameters (defaults follow the thesis):
%   'KernelSize_px'   : local fit window along the meridian (default 11)
%   'R2Threshold'     : R^2 quality threshold (default 0.70)
%   'SpeedLimits'     : [cmin cmax] m/s, valid speed range (default [0.5 30])
%   'AmpThreshold'    : reject windows whose mean phasor amplitude is below
%                       this fraction of the 99th percentile amplitude
%                       (default 0.05)
%   'ExcludeCenterRadius_mm' : reject estimates closer than this to the
%                       excitation point (default 0)
%   'GridSize'        : en face grid size [Ny Nx] or scalar (default nPosPerMer)
%   'InterpMethod'    : 'natural' (default), 'linear' or 'nearest'
%   'PhasorMethod'    : 'dtft' (default) or 'fftbin'
%   'TemporalHann'    : default true
%   'Regularize'      : also compute 1D PG-TV (default true)
%   'MuTV'            : TV weight [rad/m]. Empty (default) -> MuTVRel*median(K)
%   'MuTVRel'         : relative TV weight (default 0.1; thesis sims used
%                       mu = 250 rad/m for K ~ 2200-3700 rad/m)
%   'TVIterations', 'TVTau', 'TVSigma', 'TVTheta' : 500, 0.25, 0.25, 1.0
%   'TVBoundary'      : 'neumann' (default) or 'legacy'
%   'TV_UseR2Mask'    : restrict TV support to R^2 >= R2Threshold (default false)
%   'Verbose'         : default true
%
% Output struct R (main fields):
%   x_mm, y_mm                  en face axes
%   Elastogram_1DPG             1D PG speed map (unmasked) [m/s]
%   Elastogram_1DPG_Masked      1D PG speed map, R^2 and speed masked
%   Elastogram_1DPG_TV          1D PG-TV speed map (unmasked)
%   Elastogram_1DPG_TV_Masked   1D PG-TV speed map, R^2 and speed masked
%   SparseElastogram            nearest-pixel sparse map before interpolation
%   K_map, K_TV_map, R2_map     interpolated |k|, regularized |k|, R^2
%   Kx_map, Ky_map              outward k-vector field (radial assumption)
%   Speed_Meridian, K_Meridian, R2_Meridian, Amp_Meridian : [nPos x nMer]
%   Phasor_Meridian             [nPos x nMer] complex harmonic field
%   SupportMask, TVMask, TVInfo, PhasorInfo, Params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options ==============================

[nPos, Nt, nMer] = size(Meridian_Frames);

p = inputParser;
p.FunctionName = 'OCE_PG1D_AlongMeridians';
addParameter(p, 'KernelSize_px', 11);
addParameter(p, 'R2Threshold', 0.70);
addParameter(p, 'SpeedLimits', [0.5 30]);
addParameter(p, 'AmpThreshold', 0.05);
addParameter(p, 'ExcludeCenterRadius_mm', 0);
addParameter(p, 'GridSize', nPos);
addParameter(p, 'InterpMethod', 'natural');
addParameter(p, 'PhasorMethod', 'dtft');
addParameter(p, 'TemporalHann', true);
addParameter(p, 'Regularize', true);
addParameter(p, 'MuTV', []);
addParameter(p, 'MuTVRel', 0.1);
addParameter(p, 'TVIterations', 500);
addParameter(p, 'TVTau', 0.25);
addParameter(p, 'TVSigma', 0.25);
addParameter(p, 'TVTheta', 1.0);
addParameter(p, 'TVBoundary', 'neumann');
addParameter(p, 'TV_UseR2Mask', false);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
opts = p.Results;

if nPos ~= Geom.nPosPerMer || nMer ~= Geom.nMer
    error('Meridian_Frames is [%d x %d x %d] but Geom describes %d positions x %d meridians.', ...
        nPos, Nt, nMer, Geom.nPosPerMer, Geom.nMer);
end

omega0 = 2*pi*f_exc;
Kmin = omega0 / opts.SpeedLimits(2);
Kmax = omega0 / opts.SpeedLimits(1);

L = opts.KernelSize_px;
if mod(L,2) == 0
    L = L + 1;
    warning('KernelSize_px was even; using L = %d instead.', L);
end
halfL = floor(L/2);

if opts.Verbose
    fprintf('\nOCE_PG1D_AlongMeridians()\n');
    fprintf('Meridian_Frames: [%d pos x %d t x %d meridians]\n', nPos, Nt, nMer);
    fprintf('f_exc = %.2f Hz, fs = %.2f Hz, record = %.2f excitation periods\n', ...
        f_exc, fs, Nt*f_exc/fs);
    fprintf('Kernel = %d px = %.3f mm along the meridian\n', L, L*Geom.ds_mm);
end

%% ===================== Harmonic phasor per trace =========================

[Phasor_Meridian, PhasorInfo] = OCE_HarmonicPhasor(Meridian_Frames, 2, fs, f_exc, ...
    'Method', opts.PhasorMethod, 'TemporalHann', opts.TemporalHann);

Phasor_Meridian = reshape(Phasor_Meridian, nPos, nMer);

Amp_Meridian = abs(Phasor_Meridian);
ampRef = prctile(Amp_Meridian(isfinite(Amp_Meridian)), 99);
ampThresh = opts.AmpThreshold * ampRef;

%% ===================== 1D Phase Gradient per meridian ====================

Speed_Meridian = nan(nPos, nMer);
K_Meridian     = nan(nPos, nMer);
R2_Meridian    = nan(nPos, nMer);

s_m = Geom.s_mm * 1e-3;

for b = 1:nMer

    cplx_line = Phasor_Meridian(:,b);

    for i = (1+halfL):(nPos-halfL)

        idx = (i-halfL):(i+halfL);
        cwin = cplx_line(idx);

        if any(~isfinite(cwin)) || mean(abs(cwin)) < ampThresh
            continue;
        end

        % Phase relative to the window center reduces wraps; the constant
        % offset does not change the fitted slope.
        theta_u = unwrap(angle(cwin .* conj(cwin(halfL+1))));

        s_win = s_m(idx,b);
        s_c = s_win - mean(s_win);
        th_c = theta_u - mean(theta_u);

        k_est = (s_c.' * th_c) / (s_c.' * s_c);

        SSres = sum((th_c - k_est*s_c).^2);
        SStot = sum(th_c.^2) + eps;
        R2 = 1 - SSres/SStot;

        if abs(k_est) < eps
            continue;
        end

        K_Meridian(i,b)     = k_est;
        Speed_Meridian(i,b) = omega0 / abs(k_est);
        R2_Meridian(i,b)    = R2;
    end
end

if opts.ExcludeCenterRadius_mm > 0
    nearCenter = Geom.r_mm < opts.ExcludeCenterRadius_mm;
    K_Meridian(nearCenter)     = NaN;
    Speed_Meridian(nearCenter) = NaN;
    R2_Meridian(nearCenter)    = NaN;
end

Kmag_Meridian = abs(K_Meridian);

%% ===================== En face grid ======================================

[x_mm, y_mm, Xq_mm, Yq_mm] = OCE_MB_EnFaceGrid(Geom, opts.GridSize);

rq_mm = hypot(Xq_mm - Geom.Center_mm(1), Yq_mm - Geom.Center_mm(2));
DiskMask = rq_mm <= Geom.scanLength_mm/2;

%% ===================== Project estimates to en face ======================

valid = isfinite(Speed_Meridian) & isfinite(R2_Meridian) & Speed_Meridian > 0;

Xv = Geom.X_mm(valid);
Yv = Geom.Y_mm(valid);
Cv = Speed_Meridian(valid);
Rv = R2_Meridian(valid);
Kv = Kmag_Meridian(valid);

if opts.Verbose
    fprintf('Valid meridian estimates: %d / %d\n', nnz(valid), numel(valid));
end

if nnz(valid) < 3
    error('Fewer than 3 valid 1D phase-gradient estimates were obtained.');
end

SparseElastogram = accumulateNearest(x_mm, y_mm, Xv, Yv, Cv);
SparseK          = accumulateNearest(x_mm, y_mm, Xv, Yv, Kv);
SparseR2         = accumulateNearest(x_mm, y_mm, Xv, Yv, Rv);

%% ===================== Interpolate AFTER 1D estimation ===================

[Pu, ~, grp] = uniquetol([Xv Yv], 1e-9, 'ByRows', true, 'DataScale', 1);
avg = @(v) accumarray(grp, v, [], @mean);

F_c = scatteredInterpolant(Pu(:,1), Pu(:,2), avg(Cv), opts.InterpMethod, 'none');
Elastogram_1DPG = F_c(Xq_mm, Yq_mm);

F_c.Values = avg(Rv);
R2_map = F_c(Xq_mm, Yq_mm);

F_c.Values = avg(Kv);
K_map = F_c(Xq_mm, Yq_mm);

SupportMask = DiskMask & isfinite(K_map);

Elastogram_1DPG(~SupportMask) = NaN;
R2_map(~SupportMask) = NaN;
K_map(~SupportMask) = NaN;

% Outward k-vector field under the radial-propagation assumption that
% underlies the 1D PG estimator.
ur_x = (Xq_mm - Geom.Center_mm(1)) ./ max(rq_mm, eps);
ur_y = (Yq_mm - Geom.Center_mm(2)) ./ max(rq_mm, eps);
Kx_map = K_map .* ur_x;
Ky_map = K_map .* ur_y;

QualityMask = SupportMask & (R2_map >= opts.R2Threshold);

Elastogram_1DPG_Masked = Elastogram_1DPG;
Elastogram_1DPG_Masked(~QualityMask) = NaN;
Elastogram_1DPG_Masked(Elastogram_1DPG_Masked < opts.SpeedLimits(1) | ...
                       Elastogram_1DPG_Masked > opts.SpeedLimits(2)) = NaN;

%% ===================== 2D TV on interpolated K-map =======================

K_TV_map = nan(size(K_map));
Elastogram_1DPG_TV = nan(size(K_map));
Elastogram_1DPG_TV_Masked = nan(size(K_map));
TVMask = false(size(K_map));
TVInfo = struct();

if opts.Regularize

    TVMask = SupportMask & isfinite(K_map) & K_map >= Kmin & K_map <= Kmax;

    if opts.TV_UseR2Mask
        TVMask = TVMask & R2_map >= opts.R2Threshold;
    end

    mu = opts.MuTV;
    if isempty(mu)
        mu = opts.MuTVRel * median(K_map(TVMask), 'omitnan');
    end

    [K_TV_map, TVInfo] = OCE_TVDenoise_ChambollePock2D(K_map, mu, TVMask, ...
        'NIterations', opts.TVIterations, 'Tau', opts.TVTau, ...
        'Sigma', opts.TVSigma, 'Theta', opts.TVTheta, ...
        'Boundary', opts.TVBoundary);

    K_TV_map(K_TV_map < Kmin | K_TV_map > Kmax) = NaN;

    Elastogram_1DPG_TV = omega0 ./ K_TV_map;

    Elastogram_1DPG_TV_Masked = Elastogram_1DPG_TV;
    Elastogram_1DPG_TV_Masked(~QualityMask) = NaN;

    if opts.Verbose
        fprintf('1D PG-TV: mu = %.4g rad/m, %d iterations, %.1f %% of grid regularized\n', ...
            mu, opts.TVIterations, 100*nnz(TVMask)/numel(TVMask));
    end
end

%% ===================== Output ============================================

R = struct();
R.Method = '1D PG along meridians (+ TV on interpolated |k|)';

R.x_mm = x_mm;
R.y_mm = y_mm;

R.Elastogram_1DPG           = Elastogram_1DPG;
R.Elastogram_1DPG_Masked    = Elastogram_1DPG_Masked;
R.Elastogram_1DPG_TV        = Elastogram_1DPG_TV;
R.Elastogram_1DPG_TV_Masked = Elastogram_1DPG_TV_Masked;

R.SparseElastogram = SparseElastogram;
R.SparseK          = SparseK;
R.SparseR2         = SparseR2;

R.K_map    = K_map;
R.K_TV_map = K_TV_map;
R.R2_map   = R2_map;
R.Kx_map   = Kx_map;
R.Ky_map   = Ky_map;

R.Speed_Meridian  = Speed_Meridian;
R.K_Meridian      = K_Meridian;
R.R2_Meridian     = R2_Meridian;
R.Amp_Meridian    = Amp_Meridian;
R.Phasor_Meridian = Phasor_Meridian;

R.SupportMask = SupportMask;
R.QualityMask = QualityMask;
R.TVMask      = TVMask;
R.TVInfo      = TVInfo;
R.PhasorInfo  = PhasorInfo;

R.omega0 = omega0;
R.KernelSize_px = L;
R.Params = opts;

end

%% ======================================================================= %
%                          LOCAL FUNCTIONS
% ======================================================================== %

function M = accumulateNearest(x_mm, y_mm, Xv, Yv, V)
% Average scattered values into their nearest grid pixel ([Ny x Nx]).

[~, ix] = min(abs(Xv(:) - x_mm(:).'), [], 2);
[~, iy] = min(abs(Yv(:) - y_mm(:).'), [], 2);

S = accumarray([iy ix], V(:), [numel(y_mm) numel(x_mm)], @sum, 0);
C = accumarray([iy ix], 1,    [numel(y_mm) numel(x_mm)], @sum, 0);

M = nan(size(S));
M(C > 0) = S(C > 0) ./ C(C > 0);

end
