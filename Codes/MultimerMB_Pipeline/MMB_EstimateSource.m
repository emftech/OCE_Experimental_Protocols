function Src = MMB_EstimateSource(Frames, Geom, fs, f_an, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_EstimateSource
%
% Locates the excitation point (needle tip) in the en face plane from the
% measured wavefield, assuming an outward cylindrical wave from an unknown
% point r0 = (x0, y0):
%
%   V(r) = A(r) exp(-j (k |r - r0| + phi0))
%
% For every candidate (x0, y0, k) the coherence
%
%   C = | sum_i w_i u_i exp(+j k |r_i - r0|) | / sum_i w_i,
%   u_i = V_i/|V_i|,  w_i = |V_i|
%
% is evaluated (C = 1 for a perfect cylindrical wave), and the maximum is
% taken: coarse grid search, then a fine grid around the best point.
% The fitted k gives a global "cylindrical-fit" speed c = 2 pi f_an / k.
%
% The scan centre is NOT assumed to coincide with the needle: in the
% 06-10-2026 data the needle sits ~1-2 mm away from the centre of the
% meridian pattern, which makes the propagation non-cylindrical with
% respect to the meridians (thesis, Ch. 4).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   Frames : [nPos*nMer x Nt] narrowband traces (acquisition order, NaN ok)
%   Geom   : OCE_MB_MeridianGeometry struct
%   fs     : sampling frequency [Hz];  f_an : analysis frequency [Hz]
% Optional name-value:
%   'SearchHalfWidth_mm' : search square half width (default 4)
%   'kRange_radpermm'    : [kmin kmax] (default [0.2 4])
%   'MaxPoints'          : strongest samples used (default 1500)
%   'MinDist_mm'         : samples closer than this to r0 are ignored (0.8)
% Output struct Src: x0_mm, y0_mm, k_radpermm, c_mps, Coherence,
%   Coherence_Centered (same fit with r0 = [0 0]), c_Centered_mps,
%   Offset_mm (distance from the scan centre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'SearchHalfWidth_mm', 4);
addParameter(p, 'kRange_radpermm', [0.2 4]);
addParameter(p, 'MaxPoints', 1500);
addParameter(p, 'MinDist_mm', 0.8);
parse(p, varargin{:});
opts = p.Results;

P = OCE_HarmonicPhasor(Frames, 2, fs, f_an);
P = P(:);
X = Geom.X_mm(:);
Y = Geom.Y_mm(:);

ok = isfinite(P) & abs(P) > 0;
amp = abs(P(ok));
thr = 0.15*prctile(amp, 95);
ok(ok) = amp > thr;

idx = find(ok);
if numel(idx) > opts.MaxPoints
    [~, o] = sort(abs(P(idx)), 'descend');
    idx = idx(o(1:opts.MaxPoints));
end
U = P(idx)./abs(P(idx));
W = abs(P(idx));
Xo = X(idx);
Yo = Y(idx);

ks = linspace(opts.kRange_radpermm(1), opts.kRange_radpermm(2), 200);

score = @(x0, y0) evalScore(x0, y0, Xo, Yo, U, W, ks, opts.MinDist_mm);

% Coarse search
h = opts.SearchHalfWidth_mm;
best = [-inf 0 0 0];
for x0 = -h:0.25:h
    for y0 = -h:0.25:h
        [m, kb] = score(x0, y0);
        if m > best(1), best = [m x0 y0 kb]; end
    end
end

% Fine search
xc = best(2); yc = best(3);
for x0 = xc-0.3:0.05:xc+0.3
    for y0 = yc-0.3:0.05:yc+0.3
        [m, kb] = score(x0, y0);
        if m > best(1), best = [m x0 y0 kb]; end
    end
end

[mC, kC] = score(0, 0);

Src = struct();
Src.x0_mm = best(2);
Src.y0_mm = best(3);
Src.Offset_mm = hypot(best(2), best(3));
Src.k_radpermm = best(4);
Src.c_mps = 2*pi*f_an/(best(4)*1e3);
Src.Coherence = best(1);
Src.Coherence_Centered = mC;
Src.c_Centered_mps = 2*pi*f_an/(kC*1e3);
Src.nPoints = numel(idx);
Src.Options = opts;

end

%% ----------------------------------------------------------------------- %
function [m, kBest] = evalScore(x0, y0, Xo, Yo, U, W, ks, minD)
d = hypot(Xo - x0, Yo - y0);
use = d >= minD;
sc = abs(sum(W(use).*U(use).*exp(1j*d(use)*ks), 1)) / sum(W(use));
[m, i] = max(sc);
kBest = ks(i);
end
