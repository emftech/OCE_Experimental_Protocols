%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_Validate
%
% Offline check of the hardware SBM waveform BEFORE using it on the Twin
% OCT. Run it on any PC with MATLAB:
%   1. builds the waveform from SBM_TwinOCT_Config (here with the 12x12 mm
%      FOV written directly in volts = mm, i.e. 1 V/mm),
%   2. checks the synchronization constraints on the GLOBAL A-line clock,
%   3. plays a synthetic non-cylindrical harmonic wave (source off-centre
%      at the 06-10 needle position) through the exact timeline (galvo lag =
%      B-frame Signal Delay), records it the way the ATS9360 would, sorts it
%      with SBM_TwinOCT_Reorganize and compensates with SBM_TwinOCT_Recover,
%   4. compares raw / recovered phasors with the ground truth and saves the
%      figures in ./Validation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                                % SBM_Protocol root (Config, Waveform)
outDir = fullfile(fileparts(here), 'Validation');   % SBM_Protocol/Validation
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% 1) Waveform ------------------------------------------------------------
Cfg = SBM_TwinOCT_Config();
Cfg.horizontalMin_V = -6;  Cfg.horizontalMax_V = 6;     % "volts" = mm here
Cfg.verticalMin_V   = -6;  Cfg.verticalMax_V   = 6;
Cfg.HVC = 1;  Cfg.VVC = 1;
[W, Info] = SBM_TwinOCT_Waveform(Cfg);
P = Info.P;  Nphi = Info.Nphi;  Ts = Cfg.Ts_OCT;  nF = Cfg.nFrames;

fprintf('\n===== SBM Twin OCT waveform =====\n');
fprintf('P = %d samples, N_phi = %d, Ts_SBM = %.0f us\n', P, Nphi, Info.Ts_SBM*1e6);
fprintf('Frame body = %d, record R = %d A-lines (G = %d B-scans), C = %d records/frame\n', ...
    Info.Nbody, Info.R, Info.G, Info.C);
fprintf('Frame period F = %d = %d*P + %d samples (%.3f ms)\n', Info.F, Info.K, Nphi, Info.F*Ts*1e3);
fprintf('Acquisition window = %.1f ms, AO buffer = %d samples x %d ch (%.1f ms)\n', ...
    Info.AcqWindow_s*1e3, Info.Ntot, size(W,1), Info.Duration_s*1e3);
disp(Info.LabVIEW)

%% 2) Synchronization checks ----------------------------------------------
chk = struct();
chk.FmodP_equals_Nphi = mod(Info.F, P) == Nphi;
chk.FramePhaseStep    = all(mod(diff(Info.FrameStart), P) == Nphi);
eqT = sort(mod(Info.FramePhaseSamples, P));               % equivalent-time positions
chk.FramesCoverNPperiods = isequal(unique(mod((0:nF-1)*Nphi, P)), ...
    unique(mod(eqT - eqT(1), P))');
chk.RecordsDoNotOverlap = Info.C*Info.R <= Info.F;
rise = find(diff([0; W(3,:)'] > 2.5) == 1) - 1;            % 0-based rising edges
chk.TriggerCount = numel(rise) == Info.C*nF;
chk.TriggersAsPlanned = isequal(rise(:), Info.TrigIdx(:));
chk.WithinAORange = max(abs(W(:))) <= 10;
chk.LengthMultipleOfP = mod(Info.Ntot, P) == 0;
if size(W,1) == 4 && Cfg.ExcMode == "sine"
    % excitation phase at every frame start advances by exactly N_phi samples
    phStart = mod(Info.FrameStart, P);
    chk.ExcPhaseAtFrameStart = all(mod(diff(phStart), P) == Nphi);
    chk.ExcSineFrequency_Hz = 1/(P*Ts);
end
disp(chk)
flags = struct2cell(chk);  flags = flags(cellfun(@islogical, flags));
assert(all([flags{:}]), 'A synchronization check failed (see chk).');

%% 3) Acquisition-vector figures (thesis Fig. 49/50 style) ----------------
t_ms = (0:Info.Ntot-1)*Ts*1e3;
f1 = figure('Color','w','Position',[60 60 1100 720]);
tl = tiledlayout(f1, size(W,1), 1, 'TileSpacing','compact');
names = {'X galvo [V]','Y galvo [V]','B-trigger [V]','Excitation [V]'};
win = [Info.FrameStart(1)-2*P, Info.FrameStart(3)]*Ts*1e3;
for r = 1:size(W,1)
    nexttile(tl); plot(t_ms, W(r,:), 'LineWidth', 1); grid on;
    ylabel(names{r}); xlim(win);
end
xlabel(tl, 'Time [ms]');
title(tl, sprintf('SBM AO waveform, first two frames (N_B = %d, R = %d A-lines/trigger)', ...
    Cfg.nBscans, Info.R));
exportgraphics(f1, fullfile(outDir, 'SBM_AO_Waveform_TwoFrames.png'), 'Resolution', 150);

f2 = figure('Color','w','Position',[80 80 1100 760]);
tl2 = tiledlayout(f2, 4, 1, 'TileSpacing','compact');
for m = 1:4
    nexttile(tl2);
    i0 = Info.FrameStart(m);
    idx = (i0 - P):(i0 + Cfg.TrigDelay + 2*P);
    tt = (idx - i0)*Ts*1e6;
    yyaxis left; plot(tt, W(1,idx+1), 'b-', 'LineWidth', 1.2); ylabel('X [V]');
    yyaxis right;
    if size(W,1) == 4
        plot(tt, W(4,idx+1)/max(abs(W(4,:))), 'k-', 'LineWidth', 1.2); hold on;
    end
    stairs(tt, W(3,idx+1)/Cfg.TrigAmp_V, 'g-', 'LineWidth', 1.5);
    xline(0, ':'); grid on; ylim([-1.2 1.2]);
    title(sprintf('Frame %d: excitation phase at frame start = %d samples', m, ...
        mod(i0, P)));
end
xlabel(tl2, 'Time relative to frame start [\mus]');
exportgraphics(f2, fullfile(outDir, 'SBM_FrameWise_PhaseStepping.png'), 'Resolution', 150);

f3 = figure('Color','w','Position',[100 100 620 600]);
body = Info.FrameStart(1) + (1:Info.Nbody);
plot(W(1,body), W(2,body), '-', 'Color', [.8 .8 .8]); hold on;
isA = Info.IsAlineBody;
plot(W(1,body(~isA)), W(2,body(~isA)), '.', 'Color', [.6 0 0]);
plot(W(1,body(isA)), W(2,body(isA)), 'bo', 'MarkerFaceColor','b','MarkerSize',3);
fb = Info.FrameStart(1) + Info.Nbody + (1:Cfg.FlybackMinSamples);
plot(W(1,fb), W(2,fb), 'g-', 'LineWidth', 1.5);
axis equal; grid on; xlabel('X [V]'); ylabel('Y [V]');
legend('Raster path','Turn samples','Acquired A-lines','Inter-frame flyback','Location','southoutside');
title(sprintf('%d A-lines x %d B-scans, %d turn samples', Cfg.nAscans, Cfg.nBscans, Cfg.nTurnSamples));
exportgraphics(f3, fullfile(outDir, 'SBM_RasterPattern.png'), 'Resolution', 150);

%% 4) Synthetic acquisition through the real timeline ---------------------
c0 = 10;  f0 = Info.f_exc;  w0 = 2*pi*f0;  k0 = w0/(c0*1e3);   % rad/mm
src = [0.70 1.10];                                         % needle tip (06-10)
field = @(x,y,n) exp(-hypot(x-src(1),y-src(2))/8) .* ...
                 sin(w0*n*Ts - k0*hypot(x-src(1),y-src(2)));

Ntot = Info.TotalAlines;
acqN = zeros(1, Ntot);  ptr = 0;
for m = 1:nF
    for c = 1:Info.C
        acqN(ptr+(1:Info.R)) = Info.TrigIdx(c,m) + Cfg.AlazarLatency + (0:Info.R-1);
        ptr = ptr + Info.R;
    end
end
cmd = acqN - Cfg.TrigDelay - Cfg.AlazarLatency;            % galvo lag = trigger delay
xb = W(1, cmd+1);  yb = W(2, cmd+1);
dataVec = field(xb, yb, acqN);                             % what the OCT phase would see

[Frames, Ax] = SBM_TwinOCT_Reorganize(dataVec, Info);
[Rec, U, PhiAcq, Uraw] = SBM_TwinOCT_Recover(Frames, Info);

[XX, YY] = meshgrid(Ax.x_mm, Ax.y_mm);
rr = hypot(XX-src(1), YY-src(2));
Uideal = exp(-rr/8) .* exp(1j*(-k0*rr - pi/2));            % sin -> cos(.-pi/2)
errRec = angle(U .* conj(Uideal));
errRaw = angle(Uraw .* conj(Uideal));
fprintf('Phase error vs ideal: raw RMS = %.3f rad, recovered RMS = %.2e rad\n', ...
    rms(errRaw(:)), rms(errRec(:)));
fprintf('Amplitude error (recovered): %.2e\n', max(abs(abs(U(:)) - abs(Uideal(:)))));

% local wavenumber -> speed (simple phase gradient, away from the source)
dx = mean(diff(Ax.x_mm));  dy = mean(abs(diff(Ax.y_mm)));
kx = angle(U(:,3:end).*conj(U(:,1:end-2)))/(2*dx);
ky = angle(U(1:end-2,:).*conj(U(3:end,:)))/(2*dy);
kk = hypot(kx(2:end-1,:), ky(:,2:end-1));
far = rr(2:end-1,2:end-1) > 2;
cRec = w0./(median(kk(far))*1e3);
kxr = angle(Uraw(:,3:end).*conj(Uraw(:,1:end-2)))/(2*dx);
kyr = angle(Uraw(1:end-2,:).*conj(Uraw(3:end,:)))/(2*dy);
cRaw = w0./(median(hypot(kxr(2:end-1,:), kyr(:,2:end-1)), 'all')*1e3);
fprintf('Phase-gradient speed: recovered %.2f m/s, raw (uncompensated) %.2f m/s, true %.2f m/s\n', ...
    cRec, cRaw, c0);

f4 = figure('Color','w','Position',[120 120 1250 380]);
tl4 = tiledlayout(f4, 1, 4, 'TileSpacing','compact');
imgs = {angle(Uideal), angle(Uraw), PhiAcq, angle(U)};
ttl  = {'Ideal phase', 'Sampled (raw) phase', '\phi_{acq} (exact, from timeline)', 'Recovered phase'};
for i = 1:4
    nexttile(tl4); imagesc(Ax.x_mm, Ax.y_mm, imgs{i}); axis image; set(gca,'YDir','normal');
    clim([-pi pi]); colormap(SBM_PhaseCmap); title(ttl{i}); xlabel('x [mm]'); ylabel('y [mm]');
end
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'rad';
exportgraphics(f4, fullfile(outDir, 'SBM_Recovery_Phase.png'), 'Resolution', 150);

f5 = figure('Color','w','Position',[140 140 900 360]);
ib = ceil(Cfg.nBscans/2); ia = 3;                        % central B-scan, 3rd A-line
tEq = Ax.t_equiv_s*1e3;
tFine = linspace(0, tEq(end), 800);
gt = exp(-rr(ib,ia)/8)*sin(w0*tFine*1e-3 - k0*rr(ib,ia));
plot(tFine, gt, 'k-', 'LineWidth', 1.2); hold on;
plot(tEq, squeeze(Frames(ib,ia,:)), 'ro', 'MarkerSize', 4);
plot(tEq, squeeze(Rec(ib,ia,:)), 'b.-');
grid on; xlabel('Equivalent time [ms]'); ylabel('u (a.u.)');
legend('Real wave','Sampled (SBM frames)','Recovered','Location','best');
title(sprintf('Pixel (B-scan %d, A-line %d): \\phi_{acq} = %.2f rad', ib, ia, PhiAcq(ib,ia)));
exportgraphics(f5, fullfile(outDir, 'SBM_SinglePixel.png'), 'Resolution', 150);

save(fullfile(outDir, 'SBM_Validation_Results.mat'), 'Info', 'chk', 'cRec', 'cRaw', ...
    'errRaw', 'errRec');
fprintf('Figures saved in %s\n', outDir);
