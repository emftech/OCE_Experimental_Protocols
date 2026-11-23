function R = OCE_PG2D(Movie, x_mm, y_mm, fs, f_exc, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_PG2D
%
% 2D Phase-Gradient (2D PG) and TV-regularized 2D Phase-Gradient (2D PG-TV)
% speed estimation for EXPERIMENTAL en face OCE movies.
%
% Experimental counterpart of PhaseGradient2D_FromSims_v3.m and
% PhaseGradient2D_TV_FromSims.m (simulation study, thesis Sections 2.3 and
% 3.6):
%
%   1) Harmonic phasor at f_exc of every pixel (OCE_HarmonicPhasor).
%   2) Local plane fit of the unwrapped phase in an L x L kernel:
%
%          theta(x,y) ~ kx x + ky y + theta0
%          |k| = sqrt(kx^2 + ky^2),    c = omega0 / |k|
%
%   3) 2D PG-TV: TV regularization of the raw |k| map on the pixels with
%      valid support, R^2 >= R2Threshold and |k| within the speed limits:
%
%          K_TV = argmin_K 0.5 ||K - K_raw||^2 + mu TV(K),  c_TV = w0 / K_TV
%
% The movie is typically the interpolated multimeridian movie returned by
% OCE_MB_InterpolateEnFaceFrames, but any en face movie works (e.g. raster
% BM-mode or SBM data).
%
% Maps are [Ny x Nx] in MESHGRID convention (rows = y). Unlike the
% simulation code, maps are NOT cropped: the kernel border is NaN, so all
% outputs share the grid of the input movie and of the 1D PG outputs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   Movie : [Ny x Nx x Nt] real motion movie, NaN outside the valid support
%           (alternatively a [Ny x Nx] complex phasor, see 'InputIsPhasor')
%   x_mm  : [1 x Nx] x axis [mm]
%   y_mm  : [1 x Ny] y axis [mm]
%   fs    : temporal sampling frequency [Hz]
%   f_exc : excitation frequency [Hz]
%
% Optional name-value parameters (defaults follow the thesis):
%   'KernelSize_px' : L, plane-fit kernel L x L (default 11)
%   'R2Threshold'   : default 0.70
%   'SpeedLimits'   : [cmin cmax] m/s (default [0.5 30])
%   'AmpThreshold'  : reject kernels with mean amplitude below this fraction
%                     of the 99th percentile amplitude (default 0.05)
%   'SupportMask'   : additional logical [Ny x Nx] support (default all true)
%   'ExcludeCenterRadius_mm', 'Center_mm' : reject estimates close to the
%                     excitation point (default 0, [0 0])
%   'InputIsPhasor' : Movie already is the complex phasor (default false)
%   'PhasorMethod'  : 'dtft' (default) or 'fftbin'
%   'TemporalHann'  : default true
%   'Regularize'    : also compute 2D PG-TV (default true)
%   'MuTV', 'MuTVRel', 'TVIterations', 'TVTau', 'TVSigma', 'TVTheta',
%   'TVBoundary'    : same meaning as in OCE_PG1D_AlongMeridians
%   'Verbose'       : default true
%
% Output struct R (main fields):
%   Elastogram_2DPG, Elastogram_2DPG_Masked
%   Elastogram_2DPG_TV, Elastogram_2DPG_TV_Masked
%   Kx_map, Ky_map, K_map, K_TV_map, R2_map, Phasor
%   SupportMask, QualityMask, TVMask, TVInfo, PhasorInfo, Params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================ Parse options ==============================

Ny = size(Movie,1);
Nx = size(Movie,2);

p = inputParser;
p.FunctionName = 'OCE_PG2D';
addParameter(p, 'KernelSize_px', 11);
addParameter(p, 'R2Threshold', 0.70);
addParameter(p, 'SpeedLimits', [0.5 30]);
addParameter(p, 'AmpThreshold', 0.05);
addParameter(p, 'SupportMask', true(Ny, Nx));
addParameter(p, 'ExcludeCenterRadius_mm', 0);
addParameter(p, 'Center_mm', [0 0]);
addParameter(p, 'InputIsPhasor', false);
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
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
opts = p.Results;

if numel(x_mm) ~= Nx || numel(y_mm) ~= Ny
    error('x_mm must have %d elements and y_mm %d elements (Movie is [Ny x Nx x Nt]).', Nx, Ny);
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

x_m = x_mm(:).' * 1e-3;
y_m = y_mm(:) * 1e-3;

if opts.Verbose
    fprintf('\nOCE_PG2D()\n');
    fprintf('Grid [%d x %d], dx = %.4f mm, dy = %.4f mm\n', Ny, Nx, ...
        abs(x_mm(2)-x_mm(1)), abs(y_mm(2)-y_mm(1)));
    fprintf('Kernel = %d x %d px = %.3f x %.3f mm\n', L, L, ...
        L*abs(x_mm(2)-x_mm(1)), L*abs(y_mm(2)-y_mm(1)));
end

%% ===================== Harmonic phasor ===================================

if opts.InputIsPhasor
    Phasor = Movie;
    PhasorInfo = struct('Method', 'input', 'f_exc', f_exc, 'f_used', f_exc);
else
    [Phasor, PhasorInfo] = OCE_HarmonicPhasor(Movie, 3, fs, f_exc, ...
        'Method', opts.PhasorMethod, 'TemporalHann', opts.TemporalHann);
    if opts.Verbose
        fprintf('Phasor: %s at %.2f Hz, record = %.2f excitation periods\n', ...
            PhasorInfo.Method, PhasorInfo.f_used, PhasorInfo.NumExcitationPeriods);
    end
end

Amp = abs(Phasor);
ampRef = prctile(Amp(isfinite(Amp)), 99);
ampThresh = opts.AmpThreshold * ampRef;

SupportMask = logical(opts.SupportMask) & isfinite(Phasor);

%% ===================== Sliding-window 2D Phase Gradient ==================

% The kernel design matrix is identical for every window on a uniform grid,
% so its pseudo-inverse is computed once.
dx = (-halfL:halfL) * mean(diff(x_m));
dy = (-halfL:halfL).' * mean(diff(y_m));
[DX, DY] = meshgrid(dx, dy);
A = [DX(:), DY(:), ones(L*L,1)];
pinvA = pinv(A);

Kx_map = nan(Ny, Nx);
Ky_map = nan(Ny, Nx);
R2_map = nan(Ny, Nx);

for iy = (1+halfL):(Ny-halfL)

    rows = (iy-halfL):(iy+halfL);

    for ix = (1+halfL):(Nx-halfL)

        cols = (ix-halfL):(ix+halfL);

        % Skip windows that extend outside the valid support.
        if ~all(SupportMask(rows, cols), 'all')
            continue;
        end

        cwin = Phasor(rows, cols);

        if mean(abs(cwin), 'all') < ampThresh
            continue;
        end

        % Phase relative to the kernel center reduces wraps before the
        % separable unwrapping used in the thesis (columns, then rows).
        theta = angle(cwin .* conj(cwin(halfL+1, halfL+1)));
        theta_u = unwrap(theta, [], 1);
        theta_u = unwrap(theta_u, [], 2);

        coef = pinvA * theta_u(:);

        theta_hat = A * coef;
        SSres = sum((theta_u(:) - theta_hat).^2);
        SStot = sum((theta_u(:) - mean(theta_u(:))).^2) + eps;

        Kx_map(iy, ix) = coef(1);
        Ky_map(iy, ix) = coef(2);
        R2_map(iy, ix) = 1 - SSres/SStot;
    end
end

K_map = hypot(Kx_map, Ky_map);
K_map(K_map <= eps) = NaN;

if opts.ExcludeCenterRadius_mm > 0
    [Xq, Yq] = meshgrid(x_mm, y_mm);
    nearCenter = hypot(Xq - opts.Center_mm(1), Yq - opts.Center_mm(2)) < opts.ExcludeCenterRadius_mm;
    K_map(nearCenter) = NaN;
    R2_map(nearCenter) = NaN;
end

Elastogram_2DPG = omega0 ./ K_map;

QualityMask = isfinite(K_map) & R2_map >= opts.R2Threshold;

Elastogram_2DPG_Masked = Elastogram_2DPG;
Elastogram_2DPG_Masked(~QualityMask) = NaN;
Elastogram_2DPG_Masked(Elastogram_2DPG_Masked < opts.SpeedLimits(1) | ...
                       Elastogram_2DPG_Masked > opts.SpeedLimits(2)) = NaN;

if opts.Verbose
    fprintf('2D PG: %d estimated pixels, %d with R^2 >= %.2f\n', ...
        nnz(isfinite(K_map)), nnz(QualityMask), opts.R2Threshold);
end

%% ================= TV regularization on raw K magnitude ==================

K_TV_map = nan(Ny, Nx);
Elastogram_2DPG_TV = nan(Ny, Nx);
Elastogram_2DPG_TV_Masked = nan(Ny, Nx);
TVMask = false(Ny, Nx);
TVInfo = struct();

if opts.Regularize

    TVMask = isfinite(K_map) & K_map >= Kmin & K_map <= Kmax & ...
             R2_map >= opts.R2Threshold;

    mu = opts.MuTV;
    if isempty(mu)
        mu = opts.MuTVRel * median(K_map(TVMask), 'omitnan');
    end

    [K_TV_map, TVInfo] = OCE_TVDenoise_ChambollePock2D(K_map, mu, TVMask, ...
        'NIterations', opts.TVIterations, 'Tau', opts.TVTau, ...
        'Sigma', opts.TVSigma, 'Theta', opts.TVTheta, ...
        'Boundary', opts.TVBoundary);

    K_TV_map(K_TV_map < Kmin | K_TV_map > Kmax) = NaN;

    Elastogram_2DPG_TV = omega0 ./ K_TV_map;

    Elastogram_2DPG_TV_Masked = Elastogram_2DPG_TV;
    Elastogram_2DPG_TV_Masked(~QualityMask) = NaN;

    if opts.Verbose
        fprintf('2D PG-TV: mu = %.4g rad/m, %d iterations, %.1f %% of grid regularized\n', ...
            mu, opts.TVIterations, 100*nnz(TVMask)/numel(TVMask));
    end
end

%% ===================== Output ============================================

R = struct();
R.Method = '2D PG (+ TV on |k|)';

R.x_mm = x_mm;
R.y_mm = y_mm;

R.Elastogram_2DPG           = Elastogram_2DPG;
R.Elastogram_2DPG_Masked    = Elastogram_2DPG_Masked;
R.Elastogram_2DPG_TV        = Elastogram_2DPG_TV;
R.Elastogram_2DPG_TV_Masked = Elastogram_2DPG_TV_Masked;

R.Kx_map   = Kx_map;
R.Ky_map   = Ky_map;
R.K_map    = K_map;
R.K_TV_map = K_TV_map;
R.R2_map   = R2_map;
R.Phasor   = Phasor;
R.Amp      = Amp;

R.SupportMask = SupportMask;
R.QualityMask = QualityMask;
R.TVMask      = TVMask;
R.TVInfo      = TVInfo;
R.PhasorInfo  = PhasorInfo;

R.omega0 = omega0;
R.KernelSize_px = L;
R.Params = opts;

end
