function cl = MMB_SpeedCLim(MapCell, pct, step)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_SpeedCLim
%
% Common color range for a group of elastograms: percentiles pct (default
% [2 98]) of all masked pixels, rounded outward to multiples of step m/s
% (default 1). One shared range per figure keeps panels comparable, as in
% the thesis grids.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(pct), pct = [2 98]; end
if nargin < 3 || isempty(step), step = 1; end

if ~iscell(MapCell)
    MapCell = {MapCell};
end

v = [];
for i = 1:numel(MapCell)
    m = MapCell{i};
    if isempty(m), continue; end
    v = [v; m(isfinite(m))]; %#ok<AGROW>
end

if isempty(v)
    cl = [0 1];
    return;
end

cl = prctile(v, pct);
cl = [step*floor(cl(1)/step), step*ceil(cl(2)/step)];
if cl(2) <= cl(1)
    cl(2) = cl(1) + step;
end

end
