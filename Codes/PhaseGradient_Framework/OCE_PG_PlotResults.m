function figs = OCE_PG_PlotResults(Results, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_PG_PlotResults
%
% Figures for the output of OCE_PG_ProcessMB:
%
%   1) 2 x 2 panel of masked elastograms: 1D PG, 1D PG-TV, 2D PG, 2D PG-TV
%      (same layout as thesis Figures 33-34)
%   2) Diagnostics: harmonic wavefield Re{V}, amplitude |V|, and the R^2
%      maps of the 1D and 2D estimators
%   3) 2D PG k-vector field over |k|, with the multimeridian pattern
%
% Optional name-value parameters:
%   'CLim'         : speed color limits [m/s] (default: 5-95 percentiles)
%   'Colormap'     : default 'turbo'
%   'ShowPattern'  : overlay the sampled meridians (default true)
%   'Inc_Mask'     : logical ROI whose outline is drawn (default [])
%   'Bkg_Mask'     : logical ROI whose outline is drawn (default [])
%   'VectorStride' : quiver decimation in px (default 6)
%   'SaveDir'      : if not empty, save .fig and .png in this folder
%   'Tag'          : string appended to the figure file names (default '')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.FunctionName = 'OCE_PG_PlotResults';
addParameter(p, 'CLim', []);
addParameter(p, 'Colormap', 'turbo');
addParameter(p, 'ShowPattern', true);
addParameter(p, 'Inc_Mask', []);
addParameter(p, 'Bkg_Mask', []);
addParameter(p, 'VectorStride', 6);
addParameter(p, 'SaveDir', '');
addParameter(p, 'Tag', '');
parse(p, varargin{:});
opts = p.Results;

x = Results.x_mm;
y = Results.y_mm;
Geom = Results.Geom;

titles = struct('PG1D', '1D PG', 'PG1D_TV', '1D PG-TV', 'PG2D', '2D PG', 'PG2D_TV', '2D PG-TV');
names = fieldnames(Results.Maps);

clim_ = opts.CLim;
if isempty(clim_)
    allv = [];
    for i = 1:numel(names)
        v = Results.Maps.(names{i});
        allv = [allv; v(isfinite(v))]; %#ok<AGROW>
    end
    clim_ = prctile(allv, [5 95]);
    if ~all(isfinite(clim_)) || clim_(1) >= clim_(2)
        clim_ = [0 10];
    end
end

figs = gobjects(0);

%% ===================== 1) Elastogram panel ===============================

figs(end+1) = figure('Name', 'OCE PG Elastograms', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.08 0.08 0.62 0.8]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

order = {'PG1D', 'PG1D_TV', 'PG2D', 'PG2D_TV'};
for i = 1:numel(order)
    nexttile;
    if ~isfield(Results.Maps, order{i})
        axis off;
        continue;
    end
    E = Results.Maps.(order{i});
    v = E(isfinite(E));
    showMap(x, y, E, clim_, opts.Colormap);
    decorate(Geom, opts);
    title({titles.(order{i}), sprintf('median %.2f m/s, IQR [%.2f, %.2f]', ...
        median(v), prctile(v,25), prctile(v,75))});
    cb = colorbar;
    ylabel(cb, 'Phase speed [m/s]');
end

sgtitle(sprintf('Multimeridian MB-mode elastograms (%d meridians, f_{exc} = %.0f Hz, R^2 \\geq %.2f)', ...
    Geom.nMer, Results.f_exc, Results.Params.R2Threshold));

%% ===================== 2) Diagnostics ====================================

figs(end+1) = figure('Name', 'OCE PG Diagnostics', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.62 0.8]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

if isfield(Results, 'PG2D')
    V = Results.PG2D.Phasor;
    Vn = V / prctile(abs(V(isfinite(V))), 99);

    nexttile;
    showMap(x, y, real(Vn), [-1 1], 'parula');
    decorate(Geom, opts);
    title('Harmonic wavefield Re\{V(x,y,\omega_0)\} (interpolated)');
    colorbar;

    nexttile;
    showMap(x, y, abs(Vn), [0 1], 'parula');
    decorate(Geom, opts);
    title('Harmonic amplitude |V| (normalized)');
    colorbar;
end

if isfield(Results, 'PG1D')
    nexttile;
    showMap(x, y, Results.PG1D.R2_map, [0 1], 'parula');
    decorate(Geom, opts);
    title('1D PG fit R^2 (interpolated)');
    colorbar;
end

if isfield(Results, 'PG2D')
    nexttile;
    showMap(x, y, Results.PG2D.R2_map, [0 1], 'parula');
    decorate(Geom, opts);
    title('2D PG plane-fit R^2');
    colorbar;
end

%% ===================== 3) k-vector field =================================

if isfield(Results, 'PG2D')
    figs(end+1) = figure('Name', 'OCE PG k-vector field', 'Color', 'w');

    K = Results.PG2D.K_map;
    K(~Results.PG2D.QualityMask) = NaN;
    showMap(x, y, K, [], opts.Colormap);
    cb = colorbar;
    ylabel(cb, '|k| [rad/m]');
    hold on;

    s = opts.VectorStride;
    [Xq, Yq] = meshgrid(x, y);
    Kx = Results.PG2D.Kx_map;
    Ky = Results.PG2D.Ky_map;

    % The sign of the fitted phase slope depends on the Fourier convention;
    % vectors are oriented away from the excitation point for readability.
    outward = sign(Kx.*(Xq - Geom.Center_mm(1)) + Ky.*(Yq - Geom.Center_mm(2)));
    outward(outward == 0) = 1;
    Kx = Kx .* outward;
    Ky = Ky .* outward;
    Kx(~Results.PG2D.QualityMask) = NaN;
    Ky(~Results.PG2D.QualityMask) = NaN;

    idy = 1:s:numel(y);
    idx = 1:s:numel(x);
    quiver(Xq(idy,idx), Yq(idy,idx), Kx(idy,idx), Ky(idy,idx), 'k', 'LineWidth', 1);
    decorate(Geom, opts);
    hold off;
    title('2D PG local wave-vector field (R^2 masked)');
end

%% ===================== Save ==============================================

if ~isempty(opts.SaveDir)
    if ~exist(opts.SaveDir, 'dir')
        mkdir(opts.SaveDir);
    end
    fileNames = {'PG_Elastograms', 'PG_Diagnostics', 'PG_KVectorField'};
    for i = 1:numel(figs)
        base = fullfile(opts.SaveDir, [fileNames{i} opts.Tag]);
        savefig(figs(i), [base '.fig']);
        exportgraphics(figs(i), [base '.png'], 'Resolution', 200);
    end
end

end

%% ======================================================================= %
%                          LOCAL FUNCTIONS
% ======================================================================== %

function showMap(x, y, M, clim_, cmap)

h = imagesc(x, y, M);
set(h, 'AlphaData', isfinite(M));
set(gca, 'Color', [0.85 0.85 0.85]);
axis xy image;
colormap(gca, cmap);
if ~isempty(clim_)
    clim(clim_);
end
xlabel('X [mm]');
ylabel('Y [mm]');

end

function decorate(Geom, opts)

hold on;

if opts.ShowPattern
    for b = 1:Geom.nMer
        plot(Geom.X_mm([1 end],b), Geom.Y_mm([1 end],b), ':', ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 0.75);
    end
end

if ~isempty(opts.Inc_Mask) && any(opts.Inc_Mask(:))
    contour(xlim_axis(gca,1), xlim_axis(gca,2), double(opts.Inc_Mask), [0.5 0.5], ...
        'w', 'LineWidth', 2);
end

if ~isempty(opts.Bkg_Mask) && any(opts.Bkg_Mask(:))
    contour(xlim_axis(gca,1), xlim_axis(gca,2), double(opts.Bkg_Mask), [0.5 0.5], ...
        'k--', 'LineWidth', 1.5);
end

hold off;

end

function v = xlim_axis(ax, which)
% Axis vectors of the image displayed in ax (x for which = 1, y for 2).

img = findobj(ax, 'Type', 'image');
img = img(end);
if which == 1
    v = linspace(img.XData(1), img.XData(end), size(img.CData, 2));
else
    v = linspace(img.YData(1), img.YData(end), size(img.CData, 1));
end

end
