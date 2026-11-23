function [Movie, Interp] = OCE_MB_InterpolateEnFaceFrames(Meridian_Frames, Geom, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_MB_InterpolateEnFaceFrames
%
% Experimental counterpart of Prepare_MB_2DPhaseGradient.m (thesis Section
% 3.6.2): every sparse multimeridian en face frame is interpolated onto a
% dense Cartesian grid, so that the 2D phase-gradient estimator can compute
% both Cartesian components of the local wave vector.
%
% Differences with the simulation version:
%   - Interpolation is performed on the EXACT physical (x,y) position of
%     every MB sample instead of on their nearest reduced-grid pixel.
%   - Coincident samples (e.g. the center shared by all meridians) are
%     merged by averaging their traces.
%   - The triangulation is built once and reused for all frames.
%
% Output grid uses MESHGRID convention: Movie is [Ny x Nx x Nt].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   Meridian_Frames : [nPosPerMer x Nt x nMer]
%   Geom            : struct from OCE_MB_MeridianGeometry
%
% Optional name-value parameters:
%   'GridSize'     : scalar N or [Ny Nx] (default nPosPerMer)
%   'InterpMethod' : 'linear' (Delaunay, default), 'natural' or 'nearest'
%   'Smoothing'    : NaN-aware Gaussian smoothing of every frame
%                    (default true, as in the simulation live scripts)
%   'Sigma_px'     : Gaussian sigma in pixels (default 0.8)
%   'Verbose'      : default true
%
% Outputs:
%   Movie  : [Ny x Nx x Nt] interpolated en face movie (NaN outside hull)
%   Interp : struct with x_mm, y_mm, HullMask, SampleMask (pixels that
%            contain an original sample), sample coordinates and options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nPos, Nt, nMer] = size(Meridian_Frames);

p = inputParser;
p.FunctionName = 'OCE_MB_InterpolateEnFaceFrames';
addParameter(p, 'GridSize', nPos);
addParameter(p, 'InterpMethod', 'linear', @(s) any(strcmpi(s, {'linear','natural','nearest'})));
addParameter(p, 'Smoothing', true);
addParameter(p, 'Sigma_px', 0.8);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});
opts = p.Results;

%% ===================== Sample coordinates and traces =====================

Traces = reshape(permute(Meridian_Frames, [1 3 2]), nPos*nMer, Nt);
Xs = Geom.X_mm(:);
Ys = Geom.Y_mm(:);

% Drop traces with no valid samples, fill isolated NaNs along time.
validTrace = any(isfinite(Traces), 2);
Traces = Traces(validTrace,:);
Xs = Xs(validTrace);
Ys = Ys(validTrace);

tIdx = 1:Nt;
for i = find(any(~isfinite(Traces), 2)).'
    v = isfinite(Traces(i,:));
    if nnz(v) == 1
        Traces(i,:) = Traces(i,v);
    else
        Traces(i,~v) = interp1(tIdx(v), Traces(i,v), tIdx(~v), 'linear', 'extrap');
    end
end

% Merge coincident samples.
[Pu, ~, grp] = uniquetol([Xs Ys], 1e-9, 'ByRows', true, 'DataScale', 1);
nU = size(Pu,1);
TracesU = zeros(nU, Nt);
for t = 1:Nt
    TracesU(:,t) = accumarray(grp, Traces(:,t), [nU 1], @mean);
end

if nU < 3
    error('Not enough valid sampled positions to perform 2D interpolation.');
end

if opts.Verbose
    fprintf('\nOCE_MB_InterpolateEnFaceFrames()\n');
    fprintf('Valid sampled positions: %d (%d merged duplicates)\n', nU, numel(Xs) - nU);
    fprintf('Interpolation method: %s, smoothing: %d (sigma = %.2f px)\n', ...
        opts.InterpMethod, opts.Smoothing, opts.Sigma_px);
end

%% ===================== Grid and interpolation hull =======================

[x_mm, y_mm, Xq_mm, Yq_mm] = OCE_MB_EnFaceGrid(Geom, opts.GridSize);
[Ny, Nx] = size(Xq_mm);

F = scatteredInterpolant(Pu(:,1), Pu(:,2), ones(nU,1), 'linear', 'none');
HullMask = isfinite(F(Xq_mm, Yq_mm));

F.Method = lower(opts.InterpMethod);

%% ===================== Interpolate every en face frame ===================

Movie = nan(Ny, Nx, Nt);

for t = 1:Nt
    F.Values = TracesU(:,t);
    frame = F(Xq_mm, Yq_mm);
    frame(~HullMask) = NaN;
    Movie(:,:,t) = frame;
end

%% ===================== Optional NaN-aware Gaussian smoothing =============

if opts.Smoothing
    ksz = 2*ceil(3*opts.Sigma_px) + 1;
    weight = double(HullMask);
    weightBlur = imgaussfilt(weight, opts.Sigma_px, 'FilterSize', ksz, 'Padding', 0);
    for t = 1:Nt
        frame = Movie(:,:,t);
        frame(~HullMask) = 0;
        frameBlur = imgaussfilt(frame, opts.Sigma_px, 'FilterSize', ksz, 'Padding', 0);
        frameBlur = frameBlur ./ max(weightBlur, 1e-8);
        frameBlur(~HullMask) = NaN;
        Movie(:,:,t) = frameBlur;
    end
end

%% ===================== Output metadata ===================================

SampleMask = false(Ny, Nx);
[~, ix] = min(abs(Pu(:,1) - x_mm), [], 2);
[~, iy] = min(abs(Pu(:,2) - y_mm), [], 2);
SampleMask(sub2ind([Ny Nx], iy, ix)) = true;

Interp = struct();
Interp.x_mm = x_mm;
Interp.y_mm = y_mm;
Interp.HullMask = HullMask;
Interp.SampleMask = SampleMask;
Interp.SampleX_mm = Pu(:,1);
Interp.SampleY_mm = Pu(:,2);
Interp.Params = opts;

if opts.Verbose
    fprintf('Movie size: [%d x %d x %d], %d pixels inside hull\n', Ny, Nx, Nt, nnz(HullMask));
end

end
