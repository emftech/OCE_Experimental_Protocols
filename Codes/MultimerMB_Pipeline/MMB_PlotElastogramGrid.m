function fig = MMB_PlotElastogramGrid(MapCell, x_mm, y_mm, rowLabels, colTitles, CLim, ttl, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_PlotElastogramGrid
%
% Grid of phase-speed elastograms in the style of the thesis (Figs. 33-36):
% turbo colormap, masked pixels drawn with the darkest turbo colour,
% bold axes, row labels on the left ("8 / Meridians"), estimator titles on
% top, one shared colorbar "Phase Speed [m/s]" on the right.
%
% Inputs:
%   MapCell   : {nRows x nCols} cell of [Ny x Nx] speed maps (NaN = masked)
%   x_mm,y_mm : en face axes [mm]
%   rowLabels : {nRows x 1} cell (e.g. {'4','8','16','32'} -> "4 Meridians")
%   colTitles : {1 x nCols} cell (e.g. {'1D PG','1D PG-TV','2D PG','2D PG-TV'})
%   CLim      : [cmin cmax] m/s
%   ttl       : figure title
% Optional name-value:
%   'RowSuffix'   : text added below each row label (default 'Meridians')
%   'ShowMedian'  : write the median speed in each panel (default true)
%   'TileSize'    : panel size in pixels (default 330)
%   'Source_mm'   : needle tip [x0 y0] (or one per row, cell) drawn as '+'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'RowSuffix', 'Meridians');
addParameter(p, 'ShowMedian', true);
addParameter(p, 'TileSize', 330);
addParameter(p, 'Source_mm', []);         % [x0 y0] or {nRows x 1} of [x0 y0]: white '+' at the needle tip
parse(p, varargin{:});
opts = p.Results;

[nR, nC] = size(MapCell);
cmap = turbo(256);
bg = cmap(1,:);

fig = figure('Color', 'w', 'Position', [40 40 opts.TileSize*nC+260 opts.TileSize*nR+110]);
tl = tiledlayout(nR, nC, 'TileSpacing', 'compact', 'Padding', 'compact');

half = max(abs([x_mm(:); y_mm(:)]));
tk = -2*floor(half/2):2:2*floor(half/2);

for r = 1:nR
    for c = 1:nC
        ax = nexttile(tl);
        E = MapCell{r, c};
        if isempty(E)
            axis(ax, 'off');
            continue;
        end
        h = imagesc(ax, x_mm, y_mm, E);
        set(h, 'AlphaData', isfinite(E));
        set(ax, 'Color', bg);
        axis(ax, 'image');
        colormap(ax, cmap);
        clim(ax, CLim);
        set(ax, 'XTick', tk, 'YTick', tk, 'FontSize', 11, 'FontWeight', 'bold', 'TickDir', 'out');
        if r == 1
            title(ax, colTitles{c}, 'FontSize', 13, 'FontWeight', 'bold');
        end
        if c == 1
            if isempty(opts.RowSuffix)
                ylabel(ax, rowLabels{r}, 'FontSize', 13, 'FontWeight', 'bold');
            else
                ylabel(ax, {rowLabels{r}, opts.RowSuffix}, 'FontSize', 13, 'FontWeight', 'bold');
            end
        else
            ylabel(ax, 'Y [mm]', 'FontSize', 11, 'FontWeight', 'bold');
        end
        if r == nR
            xlabel(ax, 'X [mm]', 'FontSize', 12, 'FontWeight', 'bold');
        else
            set(ax, 'XTickLabel', []);
        end
        src = opts.Source_mm;
        if iscell(src), src = src{r}; end
        if ~isempty(src)
            hold(ax, 'on'); plot(ax, src(1), src(2), 'w+', 'MarkerSize', 10, 'LineWidth', 2); hold(ax, 'off');
        end
        if opts.ShowMedian
            v = E(isfinite(E));
            if ~isempty(v)
                text(ax, 0.03, 0.05, sprintf('%.2f m/s', median(v)), 'Units', 'normalized', ...
                    'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
                text(ax, 0.97, 0.05, sprintf('n = %d', numel(v)), 'Units', 'normalized', ...
                    'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
            end
        end
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.FontSize = 12;
cb.FontWeight = 'bold';
ylabel(cb, 'Phase Speed [m/s]', 'FontSize', 13, 'FontWeight', 'bold');
title(tl, ttl, 'FontSize', 16, 'FontWeight', 'bold');

end
