function cl = MMB_SymmetricCLim(X, pct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_SymmetricCLim
%
% Symmetric color limit for zero-mean wave fields (space-time maps, videos):
%   cl = percentile pct of |X| over finite samples (default pct = 99),
% so that the colormap is centred at zero and only ~1 % of the samples
% saturate. Returns a positive scalar (use clim([-cl cl])).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    pct = 99;
end

v = abs(double(X(isfinite(X))));
if isempty(v)
    cl = 1;
    return;
end
if numel(v) > 5e6
    v = v(randperm(numel(v), 5e6));
end
cl = prctile(v, pct);
if cl <= 0
    cl = max(v) + eps;
end

end
