function Out = SBM_TwinOCT_CheckExcitation(src, sidecarFile, Opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 18-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_CheckExcitation
%
% Checks that the excitation really is the SBM sine, from ONE captured
% record of it (optionally with the B-trigger, ao2). The excitation channel
% is the external generator output (ExcMode "gen_counter"/"gen_ao3", the
% recommended setup) or the DAQ ao3 sine (ExcMode "sine"). Always measure
% BEFORE the amplifier.
%   - is there a burst at all, how long, with the 5-period fade in/out,
%   - frequency (-> actual A-line rate = f * P), amplitude [Vpp], offset,
%     sine-fit quality (distortion / clipping / noise),
%   - with the B-trigger: number of triggers, frame period, and the
%     excitation phase at every trigger, which must advance by exactly
%     2*pi*N_phi/P per frame (25.26 deg for P = 57, N_phi = 4),
%   - generator modes: skipped cycles (the generator ignores a trigger when
%     its 1-cycle burst has not finished -> raise Cfg.Gen_DeadTime_s),
%   - comparison with what LabVIEW was given (sidecar).
%
% Sources (src):
%   'scope.csv'          scope export (time + channel columns; any header)
%   'capture.mat'        variables t [s], exc [V] (and optionally trig [V])
%   'daq'                live capture with the Data Acquisition Toolbox
%                        (wire ao3 -> ai0 and ao2 -> ai1 of an NI card with
%                        analog inputs; start it, then press Acquire in LabVIEW)
%   struct('t',..,'exc',..,'trig',..)
%
%   Out = SBM_TwinOCT_CheckExcitation('scope.csv', 'SBM_Waveform_LAST.mat')
%   Out = SBM_TwinOCT_CheckExcitation('scope.csv', '', struct('ExcCol',3,'TrigCol',2))
%   Out = SBM_TwinOCT_CheckExcitation('daq', '', struct('Device','Dev1','Duration',3))
%
% Options: TimeCol / ExcCol / TrigCol (CSV column numbers; default auto:
% first monotonic column = time, then the channel with a sine = exc and the
% one with 0/5 V pulses = trig), Probe (attenuation to undo, default 1),
% TrigLevel (2.5 V), Device/AIExc/AITrig/Rate/Duration for 'daq'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2, sidecarFile = ''; end
if nargin < 3, Opt = struct(); end
Opt = setDef(Opt, 'TimeCol', []);  Opt = setDef(Opt, 'ExcCol', []);  Opt = setDef(Opt, 'TrigCol', []);
Opt = setDef(Opt, 'Probe', 1);     Opt = setDef(Opt, 'TrigLevel', 2.5);
Opt = setDef(Opt, 'Device', 'Dev1'); Opt = setDef(Opt, 'AIExc', 'ai0'); Opt = setDef(Opt, 'AITrig', 'ai1');
Opt = setDef(Opt, 'Rate', 500e3);  Opt = setDef(Opt, 'Duration', 3);
Opt = setDef(Opt, 'MakeFigures', true);

%% Expected protocol -----------------------------------------------------------
if ~isempty(sidecarFile)
    S = load(sidecarFile);  Info = S.Info;  Wexp = S.waveform;
    fprintf('Sidecar %s: ExcMode = "%s", %d AO rows\n', sidecarFile, Info.Cfg.ExcMode, size(Wexp,1));
    needRow = ismember(lower(string(Info.Cfg.ExcMode)), ["sine", "gen_ao3", "ttl"]);
    if needRow && size(Wexp,1) < 4
        warning('SBM:NoExcRow', ['ExcMode "%s" needs a 4th AO row but the waveform LabVIEW ' ...
            'received has only %d rows.'], Info.Cfg.ExcMode, size(Wexp,1));
    end
else
    [Wexp, Info] = SBM_TwinOCT_Waveform(SBM_TwinOCT_Config());
end
Cfg = Info.Cfg;  P = Info.P;  Ts = Cfg.Ts_OCT;
mode = lower(string(Cfg.ExcMode));  isGen = ismember(mode, ["gen_counter", "gen_ao3", "gen_puff"]);
f0  = 1/Cfg.T_excitation;
VppExp = Cfg.ExcVpp_V;
stepExp = 2*pi*mod(Info.F, P)/P;                        % phase advance per frame [rad]

%% Capture -----------------------------------------------------------------------
[t, exc, trig] = getCapture(src, Opt);
exc = exc * Opt.Probe;  if ~isempty(trig), trig = trig * Opt.Probe; end
fs = 1/median(diff(t));
fprintf('Capture: %.3f s at %.3g S/s (%d points)\n', t(end)-t(1), fs, numel(t));
if fs < 20*f0
    warning('SBM:Rate', 'Sample rate %.3g S/s < 20 x f_exc: amplitude/phase will be rough.', fs);
end

%% 1. Burst detection --------------------------------------------------------------
Out = struct();  Out.f_expected_Hz = f0;  Out.Vpp_expected = VppExp;
x0  = exc - median(exc);
win = max(3, round(fs/f0));                             % one period
env = movmax(abs(x0), win);
envS = sort(env);                                      % noise = quietest 10 % of the record
noise = envS(max(1, round(0.1*numel(envS))));
thr = max(0.3*max(env), min(3*noise, 0.5*max(env)));  % works also when the sine fills the capture
act = env > thr;
Out.BurstFound = any(act) && (max(env) > 0.02);         % > 20 mV
if ~Out.BurstFound
    if isGen
        fprintf(2, ['NO SINE at the generator output (max %.1f mV).\n  Check: generator OUTPUT ' ...
            'ON, burst armed on EXTERNAL trigger, the trigger reaching EXT TRIG (counter VI ' ...
            'running / ao3 wired), the BNC to the scope.\n'], 1e3*max(abs(x0)));
        Out.Pass = false;  return
    end
    fprintf(2, ['NO SINE on the excitation channel (max %.1f mV).\n  Check: ExcMode in the ' ...
        'config / sidecar rows, the Air-puff "Transducer" = AsyncOCE (4-channel AO task) in ' ...
        'LabVIEW_Setup_Mirrors_and_Btrig.vi, the ao3 BNC, the scope channel.\n'], 1e3*max(abs(x0)));
    Out.Pass = false;  return
end
i1 = find(act, 1);  i2 = find(act, 1, 'last');
Out.BurstDuration_s = t(i2) - t(i1);
Out.BurstExpected_s = (numel(Wexp(end,:)) * Ts) * (size(Wexp,1) >= 4);
% steady part: skip the 5-period ramps + 2 periods
pad = round((Cfg.ExcRampPeriods + 2)*fs/f0);
ss = (i1+pad):(i2-pad);
if numel(ss) < 10*win, ss = i1:i2; end

%% 2. Frequency, amplitude, fit quality ---------------------------------------
xs = exc(ss);  ts = t(ss);
N  = numel(xs);  nfft = 2^nextpow2(8*N);
X  = abs(fft((xs - mean(xs)).*hann(N), nfft));
fr = (0:nfft-1)'*fs/nfft;
band = fr > 0.2*f0 & fr < 5*f0;
[~, k] = max(X.*band);
a = X(k-1); b = X(k); c = X(k+1);  dk = 0.5*(a - c)/(a - 2*b + c);
fm = (k - 1 + dk)*fs/nfft;
% least-squares sine fit at fm (refined by a small 1-D search)
fgrid = fm + linspace(-1, 1, 201)*fs/nfft;
best = inf;
for fq = fgrid
    M = [cos(2*pi*fq*ts) sin(2*pi*fq*ts) ones(N,1)];
    cf = M \ xs;  r = xs - M*cf;  e = sum(r.^2);
    if e < best, best = e;  fbest = fq;  cbest = cf;  rbest = r; end
end
Out.f_measured_Hz = fbest;
Out.f_error_ppm = (fbest/f0 - 1)*1e6;
Out.AlineRate_Hz = fbest*P;                              % the DAQ clock = A-line trigger
Out.Vpp_measured = 2*hypot(cbest(1), cbest(2));
Out.Offset_V = cbest(3);
Out.FitR2 = 1 - sum(rbest.^2)/sum((xs - mean(xs)).^2);
Out.ResidualRMS_mV = 1e3*rms(rbest);
Out.Clipping = mean(abs(xs - median(xs)) > 0.98*max(abs(xs - median(xs)))) > 0.05;
% skipped cycles: a generator that ignores triggers sits idle (at the start
% level) for whole periods. Fraction of time within +-5 % of the peak around
% the idle level, vs ~3.2 % for a pure sine (+ the burst dead time).
Apk = max(abs(xs - cbest(3)));
Out.IdleFraction = mean(abs(xs - cbest(3)) < 0.05*Apk);
if isfield(Cfg, 'Gen_DeadTime_s')
    deadFrac = Cfg.Gen_DeadTime_s*f0;                   % dead time / period
elseif isfield(Cfg, 'Gen_FreqMargin')
    deadFrac = Cfg.Gen_FreqMargin;                      % older configs
else
    deadFrac = 0.01;
end
idleExp = 2/pi*asin(0.05) + isGen*deadFrac;
Out.PeriodsChecked = round((ts(end) - ts(1))*fbest);
Out.SkippedCycles = max(0, round((Out.IdleFraction - idleExp - 0.02)*Out.PeriodsChecked));

%% 3. Phase stepping at the B-trigger -------------------------------------------
Out.TriggerChecked = ~isempty(trig);
if Out.TriggerChecked
    up = find(diff(trig > Opt.TrigLevel) == 1) + 1;
    % sub-sample edge time by linear interpolation
    te = t(up-1) + (Opt.TrigLevel - trig(up-1)) ./ (trig(up) - trig(up-1)) .* (t(up) - t(up-1));
    Out.nTriggers = numel(te);
    Out.nTriggersExpected = Info.C*Cfg.nFrames;
    fe = te(1:Info.C:end);                               % first trigger of every frame
    Out.FramePeriod_ms = 1e3*(fe(end) - fe(1))/max(numel(fe)-1, 1);
    Out.FramePeriodExpected_ms = 1e3*Info.F*Ts;
    % local phase of the excitation at every frame trigger (fit per edge on
    % +-2 periods around it, so any drift of the generator is seen too)
    phE = nan(numel(fe),1);
    for i = 1:numel(fe)
        sel = abs(t - fe(i)) < 2/fbest;                 % exactly 4 excitation periods
        if nnz(sel) < 20, continue; end
        tt = t(sel) - fe(i);
        M = [cos(2*pi*fbest*tt) sin(2*pi*fbest*tt) ones(nnz(sel),1)];
        cf = M \ exc(sel);
        phE(i) = atan2(-cf(2), cf(1));                   % phase at the trigger instant
    end
    dph = angle(exp(1j*diff(phE)));
    Out.PhaseStep_deg = rad2deg(median(dph, 'omitnan'));
    Out.PhaseStepExpected_deg = rad2deg(angle(exp(1j*stepExp)));
    Out.PhaseStepErrorStd_deg = rad2deg(std(angle(exp(1j*(dph - stepExp))), 'omitnan'));
    Out.PhaseStepOK = abs(angle(exp(1j*(deg2rad(Out.PhaseStep_deg) - stepExp)))) < deg2rad(3) ...
        && Out.PhaseStepErrorStd_deg < 3;
    Out.TriggerCountOK = Out.nTriggers == Out.nTriggersExpected;
    % A-line clock seen by the B-trigger (frame period = F A-lines) vs the one
    % implied by the sine (f * P). If they differ, the sine is NOT derived from
    % the DAQ / A-line clock (e.g. a free-running function generator).
    Out.AlineRateFromTrig_Hz = Info.F * (numel(fe)-1) / (fe(end) - fe(1));  % whole span: ~ppm precision
    Out.ClockMismatch_ppm = (Out.AlineRate_Hz/Out.AlineRateFromTrig_Hz - 1)*1e6;
    Out.LockedToAlineClock = abs(Out.ClockMismatch_ppm) < 50;
    Out.TrigTimes_s = te;  Out.FramePhase_rad = phE;
end

%% Verdict ---------------------------------------------------------------------------
Out.FrequencyOK = abs(Out.f_error_ppm) < 2000;           % A-line clock within 0.2 %
Out.AmplitudeOK = abs(Out.Vpp_measured/VppExp - 1) < 0.1;
Out.ShapeOK = Out.FitR2 > 0.98 && ~Out.Clipping && Out.SkippedCycles == 0;
flags = [Out.FrequencyOK Out.AmplitudeOK Out.ShapeOK];
if Out.TriggerChecked, flags = [flags Out.PhaseStepOK Out.TriggerCountOK Out.LockedToAlineClock]; end
Out.Pass = all(flags);

fprintf('\n=== Excitation check (%s) ===\n', mode);
if mode == "gen_counter"
    fprintf('Signal present for %.1f ms of the capture (counter trigger runs continuously)\n', ...
        1e3*Out.BurstDuration_s);
else
    fprintf('Burst: %.1f ms (AO buffer %.1f ms)\n', 1e3*Out.BurstDuration_s, 1e3*Info.Duration_s);
end
fprintf('Frequency: %.2f Hz (expected %.2f at 200 kHz, %+.0f ppm) -> f*P = %.1f Hz   [%s]\n', ...
    fbest, f0, Out.f_error_ppm, Out.AlineRate_Hz, ok(Out.FrequencyOK));
fprintf('Amplitude: %.3f Vpp (expected %.3f), offset %.3f V   [%s]\n', Out.Vpp_measured, VppExp, ...
    Out.Offset_V, ok(Out.AmplitudeOK));
fprintf('Sine fit R^2 = %.4f, residual %.1f mV rms, clipping %d, skipped cycles ~%d of %d   [%s]\n', ...
    Out.FitR2, Out.ResidualRMS_mV, Out.Clipping, Out.SkippedCycles, Out.PeriodsChecked, ok(Out.ShapeOK));
if Out.SkippedCycles > 0
    fprintf(2, ['  -> the generator ignored triggers: its 1-cycle burst lasts longer than the ' ...
        'trigger period. Shorten the generator PERIOD (raise Cfg.Gen_DeadTime_s, see ' ...
        'SBM_TwinOCT_GeneratorSettings).\n']);
end
if Out.TriggerChecked
    fprintf('Triggers: %d (expected %d), frame period %.4f ms (expected %.4f)   [%s]\n', ...
        Out.nTriggers, Out.nTriggersExpected, Out.FramePeriod_ms, Out.FramePeriodExpected_ms, ...
        ok(Out.TriggerCountOK));
    fprintf('Phase step per frame: %.2f deg (expected %.2f), scatter %.2f deg   [%s]\n', ...
        Out.PhaseStep_deg, Out.PhaseStepExpected_deg, Out.PhaseStepErrorStd_deg, ok(Out.PhaseStepOK));
    fprintf(['A-line clock: %.1f Hz from the triggers vs %.1f Hz from the sine (%+.0f ppm) ' ...
        '-> sine locked to the DAQ clock: %s\n'], Out.AlineRateFromTrig_Hz, Out.AlineRate_Hz, ...
        Out.ClockMismatch_ppm, ok(Out.LockedToAlineClock));
else
    fprintf('No trigger channel: phase stepping not checked (capture ao2 too).\n');
end
if Out.TriggerChecked && ~Out.LockedToAlineClock && isGen
    fprintf(2, ['  -> the generator is NOT being re-triggered every %d A-lines: check Burst = ON, ' ...
        'Triggered, N = 1, source EXTERNAL (not internal/continuous), and that the ÷%d ' ...
        'trigger reaches EXT TRIG.\n'], P, P);
end
if Out.Pass, fprintf('VERDICT: excitation OK.\n'); else, fprintf(2, 'VERDICT: excitation NOT OK (see above).\n'); end

%% Figures ----------------------------------------------------------------------------
if Opt.MakeFigures
    f = figure('Name','SBM excitation check','Color','w','Position',[60 60 1300 760]);
    tl = tiledlayout(f, 2, 3, 'TileSpacing','compact');
    nexttile(tl, [1 2]);
    plot(t*1e3, exc, 'k'); hold on; grid on;
    if ~isempty(trig), plot(t*1e3, trig/max(5,max(trig))*max(abs(exc)), 'g'); end
    xline(t([i1 i2])*1e3, 'r--'); xlabel('Time [ms]'); ylabel('V');
    title(sprintf('Whole capture: burst %.1f ms (red), trigger scaled (green)', 1e3*Out.BurstDuration_s));
    nexttile(tl);
    sel = abs(t - ts(round(end/2))) < 2.5/f0;
    plot(t(sel)*1e6, exc(sel), 'k.-'); hold on; grid on;
    tf = linspace(min(t(sel)), max(t(sel)), 500);
    plot(tf*1e6, cbest(3) + cbest(1)*cos(2*pi*fbest*tf) + cbest(2)*sin(2*pi*fbest*tf), 'r');
    xlabel('Time [\mus]'); ylabel('V'); legend('captured', 'sine fit');
    title(sprintf('Zoom: %.1f Hz, %.3f Vpp, R^2 = %.4f', fbest, Out.Vpp_measured, Out.FitR2));
    nexttile(tl);
    plot(fr(band)/1e3, 20*log10(X(band)/max(X)), 'k'); grid on; hold on;
    xline(f0/1e3, 'r--', sprintf('expected %.1f Hz', f0));
    for h = 2:4, xline(h*f0/1e3, ':', sprintf('%dx', h)); end
    xlabel('Frequency [kHz]'); ylabel('dB'); ylim([-90 5]); title('Spectrum (harmonics = distortion)');
    nexttile(tl);
    if Out.TriggerChecked
        m = (1:numel(Out.FramePhase_rad))';
        phU = unwrap(Out.FramePhase_rad);
        plot(m, rad2deg(phU - phU(1)), 'bo-', 'MarkerSize', 3); hold on; grid on;
        plot(m, rad2deg((m-1)*angle(exp(1j*stepExp))), 'r--');
        xlabel('Frame'); ylabel('Excitation phase at frame trigger [deg]');
        legend('measured', 'expected', 'Location','best');
        title(sprintf('Step %.2f deg (expected %.2f)', Out.PhaseStep_deg, Out.PhaseStepExpected_deg));
        nexttile(tl);
        plot(m(2:end), rad2deg(angle(exp(1j*(diff(Out.FramePhase_rad) - stepExp)))), 'k.-'); grid on;
        xlabel('Frame'); ylabel('Step error [deg]'); title('Phase-step error per frame');
    else
        text(0.05, 0.5, {'No trigger channel captured:', 'record ao2 too to check the phase stepping'});
        axis off;
    end
end
end

%% ============================================================================
function [t, exc, trig] = getCapture(src, Opt)
trig = [];
if isstruct(src)
    t = src.t(:);  exc = src.exc(:);
    if isfield(src, 'trig'), trig = src.trig(:); end
    return
end
src = char(src);
if strcmpi(src, 'daq')
    d = daq("ni");  d.Rate = Opt.Rate;
    addinput(d, Opt.Device, Opt.AIExc, "Voltage");
    addinput(d, Opt.Device, Opt.AITrig, "Voltage");
    fprintf('Recording %g s from %s/%s (exc) and %s (trig): press Acquire in LabVIEW now...\n', ...
        Opt.Duration, Opt.Device, Opt.AIExc, Opt.AITrig);
    T = read(d, seconds(Opt.Duration));
    t = seconds(T.Time - T.Time(1));  exc = T{:,1};  trig = T{:,2};
    return
end
[~, ~, e] = fileparts(src);
if strcmpi(e, '.mat')
    S = load(src);  t = S.t(:);  exc = S.exc(:);
    if isfield(S, 'trig'), trig = S.trig(:); end
    return
end
% scope CSV: numeric block after an arbitrary text header
A = readmatrix(src, 'FileType', 'text');
A = A(all(isfinite(A), 2), :);
nc = size(A, 2);
tc = Opt.TimeCol;
if isempty(tc)
    tc = find(arrayfun(@(c) all(diff(A(:,c)) > 0), 1:nc), 1);
    if isempty(tc), error('SBM:CSV', 'No monotonic time column found: set Opt.TimeCol.'); end
end
t = A(:, tc);  others = setdiff(1:nc, tc);
ec = Opt.ExcCol;  gc = Opt.TrigCol;
if isempty(ec) || (isempty(gc) && numel(others) > 1)
    % trigger = the channel that is mostly at two levels (0 / 5 V)
    isPulse = arrayfun(@(c) mean(abs(A(:,c) - median(A(:,c))) < 0.1) > 0.9 && range(A(:,c)) > 2, others);
    if isempty(ec)
        ec = others(find(~isPulse, 1));
        if isempty(ec), ec = others(1); end
    end
    if isempty(gc), gc = others(find(isPulse & others ~= ec, 1)); end
end
exc = A(:, ec);
if ~isempty(gc), trig = A(:, gc); end
fprintf('CSV: time = column %d, excitation = column %d, trigger = column %s\n', tc, ec, mat2str(gc));
end

function s = ok(b)
if b, s = 'OK'; else, s = 'FAIL'; end
end

function S = setDef(S, f, v)
if ~isfield(S, f) || isempty(S.(f)), S.(f) = v; end
end
