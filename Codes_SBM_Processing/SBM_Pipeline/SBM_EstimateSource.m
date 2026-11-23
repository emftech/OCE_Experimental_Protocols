function [Src, CohMap] = SBM_EstimateSource(U, x_mm, y_mm, f_an, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_EstimateSource
%
% Excitation point (needle tip) of an SBM acquisition, the full-field
% counterpart of MMB_EstimateSource: fit of an outward cylindrical wave
%   U(r) ~ A exp(-j (k |r - r0| + phi0))
% to the RECOVERED phasor, by maximising the coherence
%   |sum w u exp(+j k |r - r0|)| / sum w ,  u = U/|U|, w = |U|
% over r0 (coarse 0.25 mm grid, then 0.05 mm) and k.
%
% Inputs:
%   U     : [N_B x N_A] recovered phasor (SBM_HarmonicRecovery), NaN = excluded
%   x_mm, y_mm : grid axes (y ascending)
%   f_an  : analysis frequency [Hz] (for the reported speed)
% Optional name-value (MMB_EstimateSource defaults):
%   'SearchHalfWidth_mm' (4), 'kRange_radpermm' ([0.2 4]),
%   'MaxPoints' (1500), 'MinDist_mm' (0.8)
% Outputs:
%   Src    : x0_mm, y0_mm, Offset_mm, k_radpermm, c_mps, Coherence,
%            Coherence_Centered, c_Centered_mps, nPoints
%   CohMap : struct (x, y, C) of the coarse coherence map (sum it over the
%            acquisitions of a session for a joint needle position)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'SearchHalfWidth_mm', 4);
addParameter(p, 'kRange_radpermm', [0.2 4]);
addParameter(p, 'MaxPoints', 1500);
addParameter(p, 'MinDist_mm', 0.8);
parse(p, varargin{:});
opts = p.Results;

[X, Y] = meshgrid(x_mm, y_mm);
P = U(:);  X = X(:);  Y = Y(:);
ok = isfinite(P) & abs(P) > 0;
amp = abs(P(ok));
ok(ok) = amp > 0.15*prctile(amp, 95);
idx = find(ok);
if numel(idx) > opts.MaxPoints
    [~, o] = sort(abs(P(idx)), 'descend');
    idx = idx(o(1:opts.MaxPoints));
end
Uo = P(idx)./abs(P(idx));  W = abs(P(idx));  Xo = X(idx);  Yo = Y(idx);
ks = linspace(opts.kRange_radpermm(1), opts.kRange_radpermm(2), 200);

h = opts.SearchHalfWidth_mm;
gx = -h:0.25:h;  gy = -h:0.25:h;
Cm = zeros(numel(gy), numel(gx));  best = [-inf 0 0 0];
for ix = 1:numel(gx)
    for iy = 1:numel(gy)
        [mm, kb] = score(gx(ix), gy(iy));
        Cm(iy, ix) = mm;
        if mm > best(1), best = [mm gx(ix) gy(iy) kb]; end
    end
end
xc = best(2);  yc = best(3);
for x0 = xc-0.3:0.05:xc+0.3
    for y0 = yc-0.3:0.05:yc+0.3
        [mm, kb] = score(x0, y0);
        if mm > best(1), best = [mm x0 y0 kb]; end
    end
end
[mC, kC] = score(0, 0);

Src = struct();
Src.x0_mm = best(2);  Src.y0_mm = best(3);
Src.Offset_mm = hypot(best(2), best(3));
Src.k_radpermm = best(4);
Src.c_mps = 2*pi*f_an/(best(4)*1e3);
Src.Coherence = best(1);
Src.Coherence_Centered = mC;
Src.c_Centered_mps = 2*pi*f_an/(kC*1e3);
Src.nPoints = numel(idx);
Src.Options = opts;
CohMap = struct('x', gx, 'y', gy, 'C', Cm);

    function [m, kBest] = score(x0, y0)
        d = hypot(Xo - x0, Yo - y0);
        use = d >= opts.MinDist_mm;
        sc = abs(sum(W(use).*Uo(use).*exp(1j*d(use)*ks), 1)) / sum(W(use));
        [m, i] = max(sc);
        kBest = ks(i);
    end
end
