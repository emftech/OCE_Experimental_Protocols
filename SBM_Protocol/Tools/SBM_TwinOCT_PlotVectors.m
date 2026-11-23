function chk = SBM_TwinOCT_PlotVectors(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_PlotVectors
%
% Visual + numerical verification of the 4 SBM acquisition vectors BEFORE
% they go to LabVIEW (X galvo, Y galvo, B-trigger, excitation).
%
% Usage:
%   SBM_TwinOCT_PlotVectors                    % from SBM_TwinOCT_Config.m
%   SBM_TwinOCT_PlotVectors(Cfg)               % from a config struct
%   SBM_TwinOCT_PlotVectors(waveform, Info)    % what the generator returned
%   SBM_TwinOCT_PlotVectors('SBM_Waveform_LAST.mat')  % a LabVIEW sidecar:
%                                              % exactly what was sent to the DAQ
%   chk = SBM_TwinOCT_PlotVectors(..., 'Frame', 3, 'GalvoBandwidth_Hz', 1000)
%
% Figures (all time axes linked; use zoom/pan, and the slider jumps to a frame):
%   1. Timeline   : the 4 channels, frame bodies shaded, A-line records shaded
%                   on the B-trigger, acquired (valid) A-lines as dots on X/Y.
%   2. Kinematics : X/Y velocity [V/ms], acceleration [V/ms^2] and jerk
%                   [V/ms^3] from the sampled commands, and the X spectrum.
%                   Smooth transitions = no steps in velocity or acceleration.
%   3. Raster     : XY path of one frame + inter-frame gap, coloured by role.
%   4. Sync       : excitation phase at each frame start (must step by N_phi),
%                   the excitation channel re-sampled at one pixel across
%                   frames (must fall on one sine over N_P periods), and the
%                   record coverage of the valid A-lines.
% Output chk: struct of pass/fail checks plus kinematic maxima.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------------------------- Inputs ------------------------------------
args = varargin;  W = [];  Info = [];
if isempty(args) || ischar(args{1}) && ~isfile(args{1}) && ...
        any(strcmpi(args{1}, {'Frame','GalvoBandwidth_Hz'}))
    [W, Info] = SBM_TwinOCT_Waveform(SBM_TwinOCT_Config());
elseif isstruct(args{1}) && isfield(args{1}, 'nFrames')
    [W, Info] = SBM_TwinOCT_Waveform(args{1});  args(1) = [];
elseif isnumeric(args{1})
    W = args{1};  Info = args{2};  args(1:2) = [];
elseif (ischar(args{1}) || isstring(args{1})) && isfile(args{1})
    S = load(args{1});  W = S.waveform;  Info = S.Info;  args(1) = [];
end
p = inputParser;
addParameter(p, 'Frame', 1);
addParameter(p, 'GalvoBandwidth_Hz', []);
parse(p, args{:});
k0 = p.Results.Frame;

Cfg = Info.Cfg;  Ts = Cfg.Ts_OCT;  fs = 1/Ts;  P = Info.P;  nF = Cfg.nFrames;
N  = size(W,2);  t = (0:N-1)*Ts*1e3;                    % ms
hasExc = size(W,1) >= 4;
sm = Info.FrameStart;  Nb = Info.Nbody;
Dg = Cfg.TrigDelay;  lat = Cfg.AlazarLatency;

% acquisition index (0-based) of every valid pixel of every frame
jv  = Info.PixelIdx(:) - 1;                             % position in frame record
acq = sm' + jv + Dg + lat;                              % [Npix x nF]
cmd = acq - Dg - lat;                                   % beam = command delayed by Dg

%% ---------------------------- Checks ------------------------------------
chk = struct();
chk.Channels = size(W,1);
chk.FmodP_equals_Nphi  = mod(Info.F, P) == Info.Nphi;
chk.FramePhaseStep     = all(mod(diff(sm), P) == Info.Nphi);
rise = find(diff([0, W(3,:) > Cfg.TrigAmp_V/2]) == 1) - 1;
chk.TriggerCount       = numel(rise) == Info.C*nF;
chk.TriggersAsPlanned  = isequal(rise(:), Info.TrigIdx(:));
recStart = Info.TrigIdx + lat;  recEnd = recStart + Info.R - 1;  % per (c,m)
covered = true;
for m = 1:nF
    a = acq(:,m);
    in = any(a >= recStart(:,m)' & a <= recEnd(:,m)', 2);
    covered = covered && all(in);
end
chk.RecordsCoverAllPixels = covered;
chk.RecordsDoNotOverlap = all(diff(sort(recStart(:))) >= Info.R);
chk.WithinAORange      = max(abs(W(:))) <= 10;
chk.LengthMultipleOfP  = mod(N, P) == 0;
chk.NoNaN              = ~any(isnan(W(:)));

dt_ms = Ts*1e3;
v = diff(W(1:2,:), 1, 2)/dt_ms;                         % V/ms
a = diff(W(1:2,:), 2, 2)/dt_ms^2;                       % V/ms^2
j = diff(W(1:2,:), 3, 2)/dt_ms^3;                       % V/ms^3
chk.MaxVel_Vpms   = max(abs(v), [], 2)';
chk.MaxAcc_Vpms2  = max(abs(a), [], 2)';
chk.MaxJerk_Vpms3 = max(abs(j), [], 2)';
% "kink" indicator: largest single-sample change of acceleration relative
% to the peak acceleration (a C2 trajectory keeps this small, ~1/N_turn)
chk.AccStepRatio  = max(abs(diff(a, 1, 2)), [], 2)' ./ max(chk.MaxAcc_Vpms2, eps);
chk.Overshoot_V   = Info.Overshoot_V;
mode = lower(string(Cfg.ExcMode));
periodic = ismember(mode, ["sine", "gen_ao3", "gen_counter", "gen_puff"]);
if periodic
    % excitation is P-periodic on the A-line clock -> frame starts must step by N_phi
    chk.ExcPhaseStepsByNphi = all(mod(diff(mod(sm, P)), P) == Info.Nphi);
end
if mode == "gen_ao3"
    chk.GenTriggerRowPresent = hasExc;
    if hasExc
        er = find(diff([0, W(4,:) > Cfg.TTL_Amp_V/2]) == 1) - 1;
        chk.GenTrigSpacingIsP = all(diff(er) == P);
        chk.GenTrigCoversAcquisition = er(1) <= sm(1) - 5*P && er(end) >= max(recEnd(:));
        chk.GenTrigCount = numel(er);
    end
end
fprintf('\n=== SBM acquisition-vector checks ===\n');  disp(chk);
flags = struct2cell(chk);  isL = cellfun(@islogical, flags);
if all([flags{isL}]), fprintf('All pass/fail checks PASSED.\n');
else, fprintf(2, 'Some checks FAILED (see fields = 0).\n'); end

%% ---------------------------- Fig 1: timeline ---------------------------
names = {'X galvo [V]','Y galvo [V]','B-trigger [V]','Excitation [V]'};
if mode == "gen_ao3", names{4} = 'Generator trigger (ao3) [V]'; end
f1 = figure('Name','SBM vectors: timeline','Color','w','Position',[40 60 1250 800]);
tl = tiledlayout(f1, size(W,1), 1, 'TileSpacing','compact','Padding','compact');
ax = gobjects(size(W,1),1);
for r = 1:size(W,1)
    ax(r) = nexttile(tl);  hold(ax(r), 'on');
    yl = [min(W(r,:)) max(W(r,:))] + [-1 1]*0.08*max(eps, range(W(r,:)));
    for m = 1:nF                                         % frame bodies
        patch(ax(r), t(sm(m)+1)+[0 Nb Nb 0]*dt_ms, yl([1 1 2 2]), [0.92 0.95 1], ...
            'EdgeColor','none','HandleVisibility','off');
    end
    if r == 3                                           % record windows
        for q = 1:numel(recStart)
            patch(ax(r), t(recStart(q)+1)+[0 Info.R Info.R 0]*dt_ms, yl([1 1 2 2]), ...
                [1 0.93 0.8], 'EdgeColor','none','HandleVisibility','off');
        end
    end
    plot(ax(r), t, W(r,:), 'k-', 'LineWidth', 1);
    if r <= 2                                           % commanded position of acquired A-lines
        plot(ax(r), t(cmd(:)+1), W(r, cmd(:)+1), 'b.', 'MarkerSize', 5);
    end
    ylim(ax(r), yl);  ylabel(ax(r), names{r});  grid(ax(r), 'on');
end
linkaxes(ax, 'x');  xlabel(tl, 'Time [ms]');
title(tl, {sprintf('Blue = frame body | orange = A-line record (delayed by the galvo lag) | dots = acquired A-lines (commanded position)'), ...
    sprintf('F = %d samples (%.3f ms), N_\\phi = %d. Slider (bottom) jumps to a frame.', Info.F, Info.F*dt_ms, Info.Nphi)});
setFrame(ax, t, sm, Info, k0, P, dt_ms);
addSlider(f1, @(kk) setFrame(ax, t, sm, Info, kk, P, dt_ms), nF, k0);

%% ---------------------------- Fig 2: kinematics -------------------------
f2 = figure('Name','SBM vectors: galvo kinematics','Color','w','Position',[80 60 1250 820]);
tl2 = tiledlayout(f2, 4, 2, 'TileSpacing','compact','Padding','compact');
lab = {'X','Y'};  ax2 = gobjects(6,1);  q = 0;
for c = 1:2
    q = q+1; ax2(q) = nexttile(tl2, c);   plot(t, W(c,:), 'k'); ylabel([lab{c} ' [V]']); grid on;
    title(sprintf('%s command', lab{c}));
    q = q+1; ax2(q) = nexttile(tl2, c+2); plot(t(2:end), v(c,:), 'b'); ylabel('v [V/ms]'); grid on;
    q = q+1; ax2(q) = nexttile(tl2, c+4); plot(t(3:end), a(c,:), 'r'); ylabel('a [V/ms^2]'); grid on;
    title(sprintf('max |a| = %.3g V/ms^2, max acc. step / max |a| = %.3f', ...
        chk.MaxAcc_Vpms2(c), chk.AccStepRatio(c)));
end
linkaxes(ax2, 'x');  xlim(ax2(1), xlim(ax(1)));
nexttile(tl2, 7);
plot(t(4:end), j(1,:), 'm'); hold on; plot(t(4:end), j(2,:), 'c'); grid on;
ylabel('jerk [V/ms^3]'); legend('X','Y'); xlabel('Time [ms]'); xlim(xlim(ax(1)));
nexttile(tl2, 8);
seg = W(1, sm(1)+1 : sm(min(nF,5))+Info.F);             % a few frames of X
Xf = abs(fft(seg - mean(seg)));  fr = (0:numel(seg)-1)*fs/numel(seg);
half = fr <= 20e3;
plot(fr(half)/1e3, 20*log10(Xf(half)/max(Xf)), 'k'); grid on; hold on;
fLine = 1/(2*Info.S*Ts);  xline(fLine/1e3, 'b--', sprintf('line rate %.0f Hz', fLine));
if ~isempty(p.Results.GalvoBandwidth_Hz)
    xline(p.Results.GalvoBandwidth_Hz/1e3, 'r-', 'galvo bandwidth', 'LineWidth', 1.5);
end
ylim([-100 5]); xlabel('Frequency [kHz]'); ylabel('X spectrum [dB]');
title('Harmonics the X galvo must follow');
title(tl2, 'Galvo kinematics (finite differences of the AO samples). Smooth = no steps in v or a');

%% ---------------------------- Fig 3: raster -----------------------------
figure('Name','SBM vectors: raster','Color','w','Position',[120 80 760 700]);
hold on;
fr0 = sm(min(k0,nF));
gapIdx = fr0 + Nb + (1:Info.F-Nb);  gapIdx = gapIdx(gapIdx < N);
body = fr0 + (1:Nb);
isA = Info.IsAlineBody;
Nd = Info.RunInSamples;
rect = [Cfg.horizontalMin_V Cfg.verticalMin_V; Cfg.horizontalMax_V Cfg.verticalMax_V];
plot(rect([1 2 2 1 1],1), rect([1 1 2 2 1],2), 'k--');
plot(W(1,body), W(2,body), '-', 'Color', [.75 .75 .75]);
plot(W(1,body(~isA)), W(2,body(~isA)), '.', 'Color', [.85 .1 .1], 'MarkerSize', 7);
plot(W(1,body(isA)), W(2,body(isA)), 'bo', 'MarkerFaceColor','b','MarkerSize',3);
plot(W(1,gapIdx(1:Nd)), W(2,gapIdx(1:Nd)), '-', 'Color', [1 .55 0], 'LineWidth', 2);
plot(W(1,gapIdx(Nd+1:end-Nd)), W(2,gapIdx(Nd+1:end-Nd)), 'g-', 'LineWidth', 1.5);
plot(W(1,gapIdx(end-Nd+1:end)), W(2,gapIdx(end-Nd+1:end)), '-', 'Color', [.6 0 .8], 'LineWidth', 2);
axis equal; grid on; xlabel('X [V]'); ylabel('Y [V]');
legend('FOV','path','turn samples','acquired A-lines','run-out','flyback','run-in', ...
    'Location','southoutside','NumColumns',4);
title(sprintf('Frame %d: %d x %d A-lines, overshoot X %.3f V / Y %.3f V', min(k0,nF), ...
    Cfg.nAscans, Cfg.nBscans, Info.Overshoot_V(1), Info.Overshoot_V(2)));

%% ---------------------------- Fig 4: synchronization --------------------
f4 = figure('Name','SBM vectors: synchronization','Color','w','Position',[160 100 1250 420]);
tl4 = tiledlayout(f4, 1, 3, 'TileSpacing','compact');
nexttile(tl4);
stairs(1:nF, mod(sm, P), 'b-o', 'MarkerSize', 3); hold on;
plot(1:nF, mod(sm(1) + (0:nF-1)*Info.Nphi, P), 'r--');
xlabel('Frame m'); ylabel('Excitation phase at frame start [samples mod P]');
title(sprintf('Must step by N_\\phi = %d', Info.Nphi)); grid on; legend('vectors','ideal','Location','best');
nexttile(tl4);
if periodic
    [~, pix] = min(jv);                                 % first acquired pixel
    e0 = acq(pix,1);
    tq = (0:nF-1)*Info.Ts_SBM*1e3;
    tFine = linspace(0, nF*Info.Ts_SBM*1e3, 2000);
    ideal = Cfg.ExcOffset_V + Cfg.ExcVpp_V/2*sin(2*pi*(e0/P + tFine/(P*dt_ms)));
    if mode == "sine" && hasExc
        eS = W(4, acq(pix,:)+1);  lab = 'ao3 sine sampled at the pixel';
    else                                                % generator re-triggered every P A-lines
        eS = Cfg.ExcOffset_V + Cfg.ExcVpp_V/2*sin(2*pi*acq(pix,:)/P);
        lab = 'generator (1 cycle per trigger, model)';
    end
    plot(tFine, ideal, 'k-'); hold on;
    plot(tq, eS, 'ro', 'MarkerSize', 4); grid on;
    xlabel('Equivalent time (m-1) T_{SBM} [ms]'); ylabel('Excitation at pixel [V]');
    legend('sine of period P A-lines', lab, 'Location', 'southoutside');
    title(sprintf('Excitation seen by one pixel across %d frames (%d periods)', nF, Cfg.nPeriods));
else
    text(0.1, 0.5, sprintf('ExcMode = "%s": excitation not P-periodic', mode)); axis off;
end
nexttile(tl4);
cov = zeros(1,nF);
for m = 1:nF
    aa = acq(:,m);
    cov(m) = mean(any(aa >= recStart(:,m)' & aa <= recEnd(:,m)', 2))*100;
end
bar(1:nF, cov, 'FaceColor', [.3 .6 .3]); ylim([0 105]); grid on;
xlabel('Frame m'); ylabel('Valid A-lines inside a record [%]');
title(sprintf('R = %d A-lines x C = %d records per frame', Info.R, Info.C));
end

%% ==========================================================================
function setFrame(ax, t, sm, Info, k, P, dt_ms)
k = max(1, min(numel(sm), round(k)));
x0 = t(sm(k)+1) - 2*P*dt_ms;  x1 = t(sm(k)+1) + Info.F*dt_ms + P*dt_ms;
xlim(ax(1), [x0 x1]);
end

function addSlider(fig, cb, nF, k0)
uicontrol(fig, 'Style','text', 'Units','normalized', 'Position',[0.01 0.005 0.06 0.025], ...
    'String','Frame', 'BackgroundColor','w');
uicontrol(fig, 'Style','slider', 'Units','normalized', 'Position',[0.07 0.005 0.3 0.025], ...
    'Min',1, 'Max',max(nF,2), 'Value',k0, 'SliderStep',[1 5]/max(nF-1,1), ...
    'Callback', @(h,~) cb(h.Value));
end
