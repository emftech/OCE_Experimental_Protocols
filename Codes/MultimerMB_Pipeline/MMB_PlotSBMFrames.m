function fig = MMB_PlotSBMFrames(Periods_us, Ts_OCT_us, SaveBase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_PlotSBMFrames
%
% Feasible numbers of synchronized BM (SBM) frames vs excitation frequency,
% adapted from numFrames_vs_FreqMap / numFrames_vs_FreqMap_Experiments for
% many periods: with T = M * T_OCT (M integer), the feasible numbers of
% frames are the divisors of M. One distinct colour and marker per period,
% frame counts >= 8 labelled, legend outside the axes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(Ts_OCT_us), Ts_OCT_us = 5; end
if nargin < 3, SaveBase = ''; end

[~, o] = sort(1e6./Periods_us);
Periods_us = Periods_us(o);
n = numel(Periods_us);
cols = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19; ...
        0.30 0.75 0.93; 0.64 0.08 0.18; 0.20 0.20 0.20; 0.90 0.30 0.70; 0.10 0.60 0.50];
mk = 'osd^v><ph';

fig = figure('Color', 'w', 'Position', [60 60 1350 650]);
ax = axes(fig); hold(ax, 'on');
h = gobjects(n, 1); leg = cell(n, 1);
for k = 1:n
    M = Periods_us(k)/Ts_OCT_us;
    N = divisors(round(M));
    f = 1e6/Periods_us(k)*ones(size(N));
    c = cols(mod(k-1, size(cols,1)) + 1, :);
    h(k) = plot(ax, f, N, mk(mod(k-1, numel(mk)) + 1), 'MarkerSize', 8, 'LineWidth', 1.3, ...
        'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
    for v = N(N >= 8)
        text(ax, f(1) + 35, v, sprintf('%d', v), 'FontSize', 10, 'Color', c, 'FontWeight', 'bold');
    end
    leg{k} = sprintf('%g \\mus (%.0f Hz): N = %s', Periods_us(k), 1e6/Periods_us(k), strjoin(string(N), ', '));
end
hold(ax, 'off'); grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'Frequency [Hz]'); ylabel(ax, 'Number of Frames');
title(ax, 'Frequency vs Number of SBM Frames for the 9 Experimental Excitation Periods');
fmin = 1e6/max(Periods_us); fmax = 1e6/min(Periods_us);
xlim(ax, [fmin - 200, fmax + 250]);
ylim(ax, [-3, max(Periods_us)/Ts_OCT_us + 8]);
set(ax, 'FontSize', 13);
legend(h, leg, 'Location', 'eastoutside', 'FontSize', 11);

if ~isempty(SaveBase)
    MMB_SaveFigure(fig, SaveBase);
end

end
