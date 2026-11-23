function M = SBM_NeedleMask(x_mm, y_mm, Source_mm, Radius_mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 20-10-2026
%
% SBM_NeedleMask
%
% Near-field (needle tip) disk used by Main_Processing_SingleExp_SBM_EMF:
% pixels closer than Radius_mm (Cfg.NeedleRadius_mm) to the needle tip.
% The script applies it twice:
%   1. raster grid (x_mm, y_mm): the recovered phasor U is set to NaN there
%      (plus the data-driven shaft pixels) BEFORE the directional filter,
%      the square-grid interpolation and the PG estimators;
%   2. square grid (xq, yq): elastogram maps and the en face movies are set
%      to NaN there again (the interpolation and the PG kernel spread into
%      the hole), so the metrics (SBM_MapStats, whole map = ROI) exclude it.
%
%   M = SBM_NeedleMask(x_mm, y_mm, Source_mm, Radius_mm)   % [numel(y) x numel(x)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y] = meshgrid(x_mm, y_mm);
M = hypot(X - Source_mm(1), Y - Source_mm(2)) < Radius_mm;
end
