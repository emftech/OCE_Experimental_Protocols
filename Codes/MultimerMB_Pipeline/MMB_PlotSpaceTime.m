function [fig, CL] = MMB_PlotSpaceTime(STs, Titles, Xaxis, Time, MeridianLength_mm, MainTitle, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% MMB_PlotSpaceTime
%
% Space-time maps (one panel per map, stacked) of the surface particle
% velocity, meridians concatenated along x. Blue-white-red colormap
% (MMB_RedBlue, white = zero motion), symmetric limits = percentile
% CLimPct of |v| of the FIRST map (same scale in all panels, filter
% transients at both ends excluded), needle / missing columns in mid-gray,
% meridian borders dotted.
%
% Inputs:
%   STs      : {1 x nPanels} cell of [nPositions x Nt] maps (NaN allowed)
%   Titles   : {1 x nPanels} panel titles
%   Xaxis    : [1 x nPositions] mm;  Time : [1 x Nt] ms
%   MeridianLength_mm : meridian length (divider lines), [] for none
%   MainTitle : figure title
% Optional name-value: 'CLimPct' (97), 'EdgeSamples' (10), 'Colormap'
%   ('redblue' | 'turbo'), 'NaNColor' ([0.6 0.6 0.6]), 'CLim' ([] -> from CLimPct)
% Outputs: fig, CL (symmetric limit used)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'CLimPct', 97);
addParameter(p, 'EdgeSamples', 10);
addParameter(p, 'Colormap', 'redblue');
addParameter(p, 'NaNColor', [0.6 0.6 0.6]);
addParameter(p, 'CLim', []);               % fixed limit (e.g. from the cropped version)
parse(p, varargin{:});
o = p.Results;

if strcmpi(o.Colormap, 'turbo')
    cmap = turbo(256); nanC = cmap(1,:); lineC = 'w';
else
    cmap = MMB_RedBlue(255); nanC = o.NaNColor; lineC = 'k';
end
e = o.EdgeSamples;
CL = o.CLim;
if isempty(CL), CL = MMB_SymmetricCLim(STs{1}(:, e+1:end-e), o.CLimPct); end

nP = numel(STs);
fig = figure('Position', [40 40 1500 450*nP], 'Color', 'w');
tl = tiledlayout(nP, 1, 'TileSpacing', 'compact');
for i = 1:nP
    ax = nexttile(tl);
    h = imagesc(ax, Xaxis, Time, STs{i}.');
    set(h, 'AlphaData', isfinite(STs{i}.'));
    set(ax, 'Color', nanC);
    colormap(ax, cmap); clim(ax, [-CL CL]);
    if ~isempty(MeridianLength_mm)
        nMer = round(Xaxis(end)/MeridianLength_mm);
        hold(ax, 'on');
        for b = 1:nMer-1, xline(ax, MeridianLength_mm*b, ':', 'Color', lineC); end
        hold(ax, 'off');
    end
    ylabel(ax, 'time (ms)'); title(ax, Titles{i});
    set(ax, 'FontSize', 13);
end
xlabel(ax, 'x-axis (mm), meridians concatenated');
cb = colorbar(ax); cb.Layout.Tile = 'east'; ylabel(cb, 'Particle velocity [mm/s]');
title(tl, MainTitle, 'FontWeight', 'bold');

end
