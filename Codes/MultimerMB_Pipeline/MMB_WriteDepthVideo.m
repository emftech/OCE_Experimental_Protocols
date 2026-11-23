function CL = MMB_WriteDepthVideo(Vel, Bmode_IntLog, Mask, Xaxis, Zaxis, Time_ms, outFile, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
% Date of last modification: 20-10-2026 (transparent red-blue overlay)
%
% MMB_WriteDepthVideo
%
% Depth-resolved particle-velocity video of all meridians (concatenated)
% drawn on top of the B-mode.
%
% Default style ('Overlay' = 'transparent', 20-10-2026):
%   * colormap blue - white - red (MMB_RedBlue), symmetric limits +-CL;
%   * the overlay opacity grows with the motion amplitude,
%         alpha = AlphaMax * min(|v|/CL, 1)^AlphaGamma ,
%     so v = 0 is fully transparent and the B-mode shows through;
%   * the overlay is drawn only where the phase estimate is reliable:
%     inside Mask (tissue below the surface, no needle) AND where the
%     smoothed B-mode >= MinBmode_dB (low-signal pixels have random phase:
%     their |v| reached 1-2 mm/s vs ~0.15 mm/s rms for the wave);
%   * CL = CLimPct percentile (default 90) of |v| over those pixels, so
%     the wave (not the phase noise) sets the color scale;
%   * frames: medfilt2 [5 3] (PP_MBmode.m) + Gaussian sigma SmoothSigma px.
%
% 'Overlay' = 'legacy' reproduces the 07-10-2026 videos (turbo, fixed
% alpha 0.5, black outside the mask, CLimPct 99).
%
% Inputs:
%   Vel          : [nZ x nX x Nt] single, filtered depth velocity [mm/s]
%   Bmode_IntLog : [nZ x nX] B-mode [dB]
%   Mask         : [nZ x nX] logical (tissue below the surface, no needle)
%   Xaxis, Zaxis, Time_ms : axes
%   outFile      : .mp4 path ('' -> draw only, e.g. Live Editor animation)
% Optional name-value:
%   'Overlay' ('transparent'), 'FrameIdx' (10:Nt), 'FrameRate' (20),
%   'Quality' (50), 'CLimPct' (90), 'CLim' ([] -> from CLimPct),
%   'MinBmode_dB' (90), 'AlphaMax' (0.95), 'AlphaGamma' (0.6),
%   'SmoothSigma' (1.2), 'Alpha' (0.5, legacy), 'IThresh' ([70 130]),
%   'MeridianLength_mm' ([] -> no divider lines), 'MainTitle' (''),
%   'Visible' ('off'), 'FigWidth' ([] -> from nX), 'XLabel',
%   'VLines' (x of dashed marker lines), 'VLineLabel'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nt = size(Vel, 3);
p = inputParser;
addParameter(p, 'Overlay', 'transparent');
addParameter(p, 'FrameIdx', 10:Nt);
addParameter(p, 'FrameRate', 20);
addParameter(p, 'Quality', 50);
addParameter(p, 'CLimPct', []);
addParameter(p, 'CLim', []);
addParameter(p, 'MinBmode_dB', 90);
addParameter(p, 'AlphaMax', 0.95);
addParameter(p, 'AlphaGamma', 0.6);
addParameter(p, 'SmoothSigma', 1.2);
addParameter(p, 'Alpha', 0.5);
addParameter(p, 'IThresh', [70 130]);
addParameter(p, 'MeridianLength_mm', []);
addParameter(p, 'MainTitle', '');
addParameter(p, 'Visible', 'off');
addParameter(p, 'FigWidth', []);          % [] -> from the number of columns
addParameter(p, 'XLabel', 'x-axis (mm), meridians concatenated');
addParameter(p, 'VLines', []);            % x positions of dashed marker lines (e.g. needle tip)
addParameter(p, 'VLineLabel', '');
parse(p, varargin{:});
opts = p.Results;
legacy = strcmpi(opts.Overlay, 'legacy');
if isempty(opts.CLimPct)
    opts.CLimPct = 90 + 9*legacy;
end

%% Overlay support and background
if legacy
    cmap = turbo(255);
    M = Mask;
    I = uint8(255*mat2gray(Bmode_IntLog, opts.IThresh));
    bgImg = ind2rgb(I, gray(255)).*M;
    bg = (1-opts.Alpha).*bgImg.*M + bgImg.*(1-M);
else
    cmap = MMB_RedBlue(255);
    M = Mask & medfilt2(Bmode_IntLog, [5 5], 'symmetric') >= opts.MinBmode_dB;
    bg = repmat(mat2gray(Bmode_IntLog, opts.IThresh), 1, 1, 3);
end

%% Color limit from the reliable pixels of a subset of frames
CL = opts.CLim;
if isempty(CL)
    sub = opts.FrameIdx(1:5:end);
    smp = Vel(:, :, sub);
    smp = smp(repmat(M, 1, 1, numel(sub)));
    CL = MMB_SymmetricCLim(smp, opts.CLimPct);
    clear smp
end

writeIt = ~isempty(outFile);
if writeIt
    vw = VideoWriter(outFile, 'MPEG-4');
    vw.FrameRate = opts.FrameRate;
    vw.Quality = opts.Quality;
    open(vw);
end

nX = size(Vel, 2);
W = opts.FigWidth;
if isempty(W), W = min(1800, 360 + 1.1*nX); end
fig = figure('Color', 'w', 'Position', [40 40 W 620], 'Visible', opts.Visible);
ax = axes(fig);
hIm = imagesc(ax, Xaxis, Zaxis, bg);
hold(ax, 'on');
if ~isempty(opts.MeridianLength_mm)
    nMer = round(Xaxis(end)/opts.MeridianLength_mm);
    for b = 1:nMer-1
        xline(ax, opts.MeridianLength_mm*b, ':', 'Color', [0.7 0.7 0.7]);
    end
end
for k = 1:numel(opts.VLines)
    xline(ax, opts.VLines(k), '--', opts.VLineLabel, 'Color', [1 0.85 0.2], 'LineWidth', 1.5, ...
        'LabelVerticalAlignment', 'top', 'LabelOrientation', 'horizontal', 'FontSize', 12);
end
hold(ax, 'off');
colormap(ax, cmap);
clim(ax, [-CL CL]);
cb = colorbar(ax);
ylabel(cb, 'Particle velocity [mm/s]');
ylabel(ax, 'z-axis (mm)');
xlabel(ax, opts.XLabel);
set(ax, 'FontSize', 13);

nC = size(cmap, 1);
for j = opts.FrameIdx
    Fr = medfilt2(double(Vel(:,:,j)), [5 3], 'symmetric');
    if legacy
        fg = ind2rgb(uint8(255*mat2gray(Fr, [-CL CL])), cmap);
        frame = opts.Alpha.*fg.*M + bg;
    else
        if opts.SmoothSigma > 0
            Fr = imgaussfilt(Fr, opts.SmoothSigma);
        end
        u = max(min(Fr/CL, 1), -1);
        fg = ind2rgb(uint8(round((u + 1)/2*(nC - 1))), cmap);
        a = opts.AlphaMax*abs(u).^opts.AlphaGamma.*M;
        frame = bg.*(1 - a) + fg.*a;
    end
    set(hIm, 'CData', frame);
    title(ax, sprintf('%s  Motion frames - filtered (t = %.2f ms)', opts.MainTitle, Time_ms(j)));
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
