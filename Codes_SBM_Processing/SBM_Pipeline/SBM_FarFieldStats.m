function [StatsFar, Ring] = SBM_FarFieldStats(Maps, x_mm, y_mm, Source_mm, R0_mm, Margin_mm, RingWidth_mm, c_ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 20-10-2026
%
% SBM_FarFieldStats
%
% Far-field metrics and ring profile of the SBM elastograms. The first
% ~1 mm outside the blanked needle disk is biased (the square-grid
% interpolation and the 11-px PG kernel straddle the hole; the MB chain on
% simulations: 17.2 -> 12.6 -> 10.85 m/s in 0.5 mm rings from 4.25 mm), so
%   StatsFar : SBM_MapStats over r >= R0_mm + Margin_mm (far-field ROI,
%              same definition as the simulations, 20-10-2026)
%   Ring     : table of the median speed in rings of RingWidth_mm from R0_mm
%              outwards (columns r_in_mm, r_out_mm, <method>_Median, <method>_N)
%              to check where the bias ends.
%
%   [StatsFar, Ring] = SBM_FarFieldStats(Res.Maps, xq, yq, Source_mm, 4.25, 1, 0.5, 10.54)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 8, c_ref = NaN; end
[X, Y] = meshgrid(x_mm, y_mm);
r = hypot(X - Source_mm(1), Y - Source_mm(2));
names = fieldnames(Maps);
Far = Maps;
for i = 1:numel(names)
    m = Far.(names{i});  m(r < R0_mm + Margin_mm) = NaN;  Far.(names{i}) = m;
end
StatsFar = SBM_MapStats(Far, c_ref);

edges = R0_mm:RingWidth_mm:max(r(:));
nR = max(numel(edges) - 1, 0);
Ring = table(edges(1:nR)', edges(2:nR+1)', 'VariableNames', {'r_in_mm', 'r_out_mm'});
for i = 1:numel(names)
    med = nan(nR, 1);  N = zeros(nR, 1);
    v = Maps.(names{i});
    for k = 1:nR
        s = v(r >= edges(k) & r < edges(k+1));  s = s(isfinite(s));
        N(k) = numel(s);  if N(k) > 0, med(k) = median(s); end
    end
    Ring.([names{i} '_Median']) = med;  Ring.([names{i} '_N']) = N;
end
end
