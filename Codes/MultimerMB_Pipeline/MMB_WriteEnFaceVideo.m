function CLused = MMB_WriteEnFaceVideo(Movies, Labels, x_mm, y_mm, Time_ms, Geom, outFile, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_WriteEnFaceVideo
%
% En face particle-velocity video of one or more interpolated movies placed
% side by side (e.g. without / with directional filtering). Since
% 20-10-2026: blue-white-red colormap (MMB_RedBlue, white = zero motion),
% symmetric color limits (percentile CLimPct of |v|), unsampled pixels and
% the needle disk in mid-gray (NaNColor), meridian lines dotted.
% 'Colormap','turbo' restores the thesis style (NaN in the darkest turbo).
%
% Inputs:
%   Movies  : {1 x nPanels} cell of [Ny x Nx x Nt] movies
%   Labels  : {1 x nPanels} panel titles
%   x_mm, y_mm, Time_ms : axes
%   Geom    : meridian geometry (dotted meridian lines), [] to skip
%   outFile : .mp4 path. If empty, frames are only drawn (Live Editor
%             animation) and nothing is written.
% Optional name-value:
%   'FrameIdx'  : frames to render (default 10:Nt)
%   'FrameRate' : 20, 'Quality' : 60, 'CLimPct' : 97, 'MainTitle' : ''
%   'Colormap'  : 'redblue' (default) or 'turbo';  'NaNColor' : [0.6 0.6 0.6]
%   'Units'     : colorbar label (default 'Particle velocity [mm/s]')
%   'CLim'      : fixed symmetric limit (e.g. the one of the cropped version,
%                 so that the near field saturates in the '_Full' video)
% Output: CLused, the symmetric color limit actually used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nt = size(Movies{1}, 3);
p = inputParser;
addParameter(p, 'FrameIdx', 10:Nt);
addParameter(p, 'FrameRate', 20);
addParameter(p, 'Quality', 60);
addParameter(p, 'CLimPct', 97);
addParameter(p, 'Colormap', 'redblue');      % 'redblue' or 'turbo'
addParameter(p, 'NaNColor', [0.6 0.6 0.6]);
addParameter(p, 'MainTitle', '');
addParameter(p, 'Units', 'Particle velocity [mm/s]');
addParameter(p, 'Visible', 'off');
addParameter(p, 'CLim', []);
parse(p, varargin{:});
opts = p.Results;

nP = numel(Movies);
if strcmpi(opts.Colormap, 'turbo')
    cmap = turbo(256); nanC = cmap(1,:); lineC = [1 1 1 0.35];
else
    cmap = MMB_RedBlue(255); nanC = opts.NaNColor; lineC = [0.25 0.25 0.25 0.4];
end
CL = zeros(1, nP);
for k = 1:nP
    CL(k) = MMB_SymmetricCLim(Movies{k}(:,:,opts.FrameIdx), opts.CLimPct);
end
CL(:) = max(CL);            % same scale in all panels -> amplitudes comparable
if ~isempty(opts.CLim), CL(:) = opts.CLim; end
CLused = CL(1);

writeIt = ~isempty(outFile);
if writeIt
    vw = VideoWriter(outFile, 'MPEG-4');
    vw.FrameRate = opts.FrameRate;
    vw.Quality = opts.Quality;
    open(vw);
end

fig = figure('Color', 'w', 'Position', [40 40 520*nP+80 560], 'Visible', opts.Visible);
tl = tiledlayout(fig, 1, nP, 'TileSpacing', 'compact', 'Padding', 'compact');
hIm = gobjects(1, nP);
for k = 1:nP
    ax = nexttile(tl);
    hIm(k) = imagesc(ax, x_mm, y_mm, Movies{k}(:,:,opts.FrameIdx(1)));
    set(ax, 'Color', nanC);
    axis(ax, 'image');
    colormap(ax, cmap);
    clim(ax, [-CL(k) CL(k)]);
    hold(ax, 'on');
    if ~isempty(Geom)
        for b = 1:Geom.nMer
            plot(ax, Geom.X_mm([1 end], b), Geom.Y_mm([1 end], b), ':', 'Color', lineC);
        end
    end
    hold(ax, 'off');
    xlabel(ax, 'X [mm]'); ylabel(ax, 'Y [mm]');
    title(ax, Labels{k});
    set(ax, 'FontSize', 13, 'FontWeight', 'bold');
end
cb = colorbar(ax);
cb.Layout.Tile = 'east';
ylabel(cb, opts.Units);
ht = title(tl, '', 'FontSize', 14, 'FontWeight', 'bold');

for j = opts.FrameIdx
    for k = 1:nP
        F = Movies{k}(:,:,j);
        set(hIm(k), 'CData', F, 'AlphaData', isfinite(F));
    end
    ht.String = sprintf('%s   t = %.2f ms', opts.MainTitle, Time_ms(j));
    drawnow;
    if writeIt
        writeVideo(vw, getframe(fig));
    end
end

if writeIt
    close(vw);
    close(fig);
end

end
