function [waveform, Info] = SBM_TwinOCT_Waveform(Cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_Waveform
%
% Hardware version of AcqVectors_SBM / AcqVectors_SBM_v2 for the Twin OCT.
% Builds the analog-output waveform that LabVIEW writes to the NI-DAQ
% (one column per A-line clock tick, Ts_OCT = 5 us):
%
%   waveform(1,:)  X galvo   [V]
%   waveform(2,:)  Y galvo   [V]
%   waveform(3,:)  B-trigger [V]  (AUX trigger-enable of the ATS9360)
%   waveform(4,:)  Excitation [V] (only if Cfg.ExcMode ~= "none")
%
% Same raster as the thesis (serpentine, y from top to bottom, Hermite
% turns with tangent 0.2), but the timing is built on the GLOBAL A-line
% clock n so that the synchronization condition holds exactly in hardware:
%
%   frame start      s_m = s_1 + (m-1) F,   F = K P + N_phi
%   => s_m mod P = s_1 + (m-1) N_phi  (mod P)          (thesis Eq. 43/54)
%
% i.e. the inter-frame segment (flyback + hold) spans an integer number of
% excitation periods plus exactly one phase step, so every frame samples
% the harmonic field N_phi samples later than the previous one.
% Excitation "sine":  e(n) = Voff + Vpp/2 sin(2 pi n / P)  (global n)
% Excitation "ttl" :  one pulse per frame at e_m = s_m - Lead - (m-1) N_phi,
%                     so e_{m+1} - e_m = K P (bursts in phase with a
%                     continuous sine of period P).
%
% Info holds everything needed to sort the recorded A-lines back into
% frames and to compute the exact accumulated acquisition phase.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(Cfg)
    Cfg = SBM_TwinOCT_Config();
end

%% ======================= Synchronization constraints ====================
Ts   = Cfg.Ts_OCT;
Praw = Cfg.T_excitation / Ts;
if abs(Praw - round(Praw)) > 1e-9
    error('SBM:Period', ['T_excitation/Ts_OCT = %.4f is not an integer. ' ...
        'Choose T_exc as a multiple of %.1f us.'], Praw, Ts*1e6);
end
P  = round(Praw);
nA = Cfg.nAscans;  nB = Cfg.nBscans;  nT = Cfg.nTurnSamples;
nF = Cfg.nFrames;  NP = Cfg.nPeriods;

NphiRaw = NP * P / nF;                                  % Eq. (46)
if abs(NphiRaw - round(NphiRaw)) > 1e-9
    error('SBM:PhaseStep', ['N_phi = N_P*P/N_F = %d*%d/%d = %.4f is not an ' ...
        'integer (Eq. 47). Feasible N_F for N_P = %d: %s'], NP, P, nF, NphiRaw, ...
        NP, mat2str(feasibleFrames(P, NP)));
end
Nphi = round(NphiRaw);
if Nphi >= P/2
    error('SBM:Nyquist', 'N_phi = %d >= P/2: the frames alias f_exc.', Nphi);
end

S = nA + nT;                                            % one B-scan step
if mod(S, P) ~= 0
    warning('SBM:BscanPhase', ['N_A + N_turn = %d is not a multiple of P = %d: ' ...
        'B-scans start at different excitation phases (still recovered ' ...
        'exactly by the accumulated phase).'], S, P);
end

%% ======================= Single-frame raster (normalized) ===============
Nbody  = nB*S - nT;
xn     = nan(Nbody,1);  yn = nan(Nbody,1);
isAl   = false(Nbody,1);
rowOf  = zeros(Nbody,1); colOf = zeros(Nbody,1);
yrow   = linspace(1, -1, nB);
for b = 1:nB
    i0  = (b-1)*S;                                      % 0-based start
    idx = i0 + (1:nA);
    if mod(b,2) == 1
        xr = linspace(-1, 1, nA);  cols = 1:nA;
    else
        xr = linspace(1, -1, nA);  cols = nA:-1:1;
    end
    xn(idx) = xr;  yn(idx) = yrow(b);
    isAl(idx) = true;  rowOf(idx) = b;  colOf(idx) = cols;
end
% B-scan turns: quintic Hermite segments that match position, velocity and
% acceleration of the linear A-line ramps on both sides (C2 continuous),
% so the galvo never sees a velocity or acceleration step. The overshoot
% beyond the FOV is the physical minimum for reversing at the scan speed
% within N_turn samples (0.3125*v*(N_turn+1) in normalized units).
vs = 2/(nA-1);                                          % ramp speed [norm./sample]
if nT == 0
    warning('SBM:NoTurn', 'N_turn = 0: the X galvo reverses instantaneously.');
end
% Long turns (N_turn >= 2*t1, e.g. the frequency sweep with S = 2P): a single
% reversal over the whole turn would overshoot proportionally to N_turn, so
% the galvo instead STOPS (quintic, t1 samples), DWELLS while Y steps, and
% GOES again (quintic, t1 samples): the overshoot stays at vs*t1/2, the
% value of the canonical 37 + 20 turn, whatever N_turn is.
Nrev = 21;
if isfield(Cfg, 'TurnReversalSamples') && ~isempty(Cfg.TurnReversalSamples)
    Nrev = Cfg.TurnReversalSamples;
end
t1 = ceil(Nrev/2);
for b = 1:(nB-1)
    tidx = (b-1)*S + nA + (1:nT);
    if isempty(tidx), continue; end
    d  = 1 - 2*(mod(b,2)==0);                           % direction of row b
    xe = xn((b-1)*S+nA);                                % row b ends where row b+1 starts
    if nT < 2*t1
        N1 = nT + 1;  s = (1:nT)'/N1;
        xn(tidx) = quintic(s, xe, d*vs*N1, xe, -d*vs*N1);
    else
        xs  = xe + d*vs*t1/2;                           % stop point
        dec = quintic((1:t1)'/t1, xe, d*vs*t1, xs, 0);
        acc = quintic((0:t1-1)'/t1, xs, 0, xe, -d*vs*t1);
        xn(tidx) = [dec; repmat(xs, nT - 2*t1, 1); acc];
    end
    yn(tidx) = quintic((1:nT)'/(nT+1), yrow(b), 0, yrow(b+1), 0);
end

Vx = [Cfg.horizontalMin_V Cfg.horizontalMax_V];
Vy = [Cfg.verticalMin_V   Cfg.verticalMax_V];
toVx = @(u) (u+1)/2*(Vx(2)-Vx(1)) + Vx(1);
toVy = @(u) (u+1)/2*(Vy(2)-Vy(1)) + Vy(1);

%% ======================= Records (Alazar) ===============================
G = Cfg.BscansPerRecord;
if isempty(G)
    G = max(1, min(nB, floor(Cfg.MaxAlinesPerBuffer / S)));
end
R = G*S;                                                % A-lines per trigger
C = ceil(nB / G);                                       % triggers per frame
if R > Cfg.MaxAlinesPerBuffer
    error('SBM:Buffer', 'Record of %d A-lines exceeds Max A-lines per buffer (%d).', ...
        R, Cfg.MaxAlinesPerBuffer);
end
Dg  = Cfg.TrigDelay;  lat = Cfg.AlazarLatency;

%% ======================= Frame period ===================================
% Inter-frame gap = run-out (decelerate from scan speed, N_d samples)
%                 + flyback (rest-to-rest, >= FlybackMinSamples)
%                 + run-in (accelerate to scan speed, N_d samples)
Lfb = Cfg.FlybackMinSamples;
Nd  = Cfg.RunInSamples;
if isempty(Nd), Nd = max(3, ceil((min(nT, Nrev-1)+1)/sqrt(2))); end  % same peak jerk as a (canonical) turn
need  = max(C*R, Nbody + 2*Nd + Lfb);                   % records may run into the gap
K     = ceil(need / P);
mode  = lower(string(Cfg.ExcMode));
Lead  = 0;  BurstCycles = 0;
if mode == "ttl"
    Lead = Cfg.TTL_LeadPeriods * P;
    BurstCycles = ceil((Lead + (nF-1)*Nphi + C*R + Dg + lat) / P) + 1;
    K = max(K, BurstCycles + Cfg.TTL_GapPeriods);
end
F = K*P + Nphi;                                         % F mod P = N_phi

pre = Cfg.PreamblePeriods * P;
if mode == "ttl"
    pre = max(pre, ceil((Lead + Cfg.TTL_Width)/P)*P);
end
s1 = pre;                                               % 0-based index of frame 1
sm = s1 + (0:nF-1)'*F;

lastAcq = sm(end) + C*R + Dg + lat;                     % exclusive end
postLen = max(Lfb, Cfg.PostPeriods*P);
if mode == "sine"                                       % room for the fade-out
    postLen = max(postLen, (Cfg.ExcRampPeriods + 2 + Cfg.PostPeriods)*P);
end
Ntot    = ceil((lastAcq + postLen) / P) * P;            % multiple of P (regeneration-safe)

%% ======================= Galvo channels =================================
% Built in normalized units (C2 everywhere: every junction matches position,
% velocity and acceleration), mapped to volts at the end.
dS = sign(xn(2) - xn(1));  dE = sign(xn(Nbody) - xn(Nbody-1));
xRun  = xn(1) - dS*vs*Nd/2;                             % run-in start (at rest)
xStop = xn(Nbody) + dE*vs*Nd/2;                         % run-out end (at rest)
sIn   = (0:Nd-1)'/Nd;                                   % last run-in sample precedes body
runIn = quintic(sIn, xRun, 0, xn(1), dS*vs*Nd);
sOut  = (1:Nd)'/Nd;
runOut = quintic(sOut, xn(Nbody), dE*vs*Nd, xStop, 0);
gap  = F - Nbody;
Lfly = gap - 2*Nd;
fly  = quintic((1:Lfly)'/Lfly, xStop, 0, xRun, 0);      % rest-to-rest (minimum jerk)
yGap = quintic((1:gap)'/(gap+1), yn(Nbody), 0, yn(1), 0);

xs = nan(Ntot,1);  ys = nan(Ntot,1);
nMove = min(s1 - Nd - 1, max(Lfb, 2*P));                % park (0,0) -> run-in start
xs(1:nMove) = quintic((1:nMove)'/nMove, 0, 0, xRun, 0);
ys(1:nMove) = quintic((1:nMove)'/nMove, 0, 0, yn(1), 0);
xs(nMove+1:s1-Nd) = xRun;  ys(nMove+1:s1) = yn(1);
xs(s1-Nd+(1:Nd)) = runIn;
for m = 1:nF
    i0 = sm(m);
    xs(i0+(1:Nbody)) = xn;  ys(i0+(1:Nbody)) = yn;
    if m < nF
        xs(i0+Nbody+(1:gap)) = [runOut; fly; runIn];
        ys(i0+Nbody+(1:gap)) = yGap;
    else                                                % run-out, then back to park
        nEnd = Ntot - (i0+Nbody);
        nBack = min(nEnd - Nd, max(Lfb, 2*P));
        xs(i0+Nbody+(1:Nd)) = runOut;
        xs(i0+Nbody+Nd+(1:nBack)) = quintic((1:nBack)'/nBack, xStop, 0, 0, 0);
        ys(i0+Nbody+(1:Nd+nBack)) = quintic((1:Nd+nBack)'/(Nd+nBack), yn(Nbody), 0, 0, 0);
        xs(i0+Nbody+Nd+nBack+1:end) = 0;  ys(i0+Nbody+Nd+nBack+1:end) = 0;
    end
end
X = toVx(xs);  Y = toVy(ys);
overshoot = [max(abs(xs)) max(abs(ys))] - 1;            % beyond the FOV (normalized)

%% ======================= B-trigger channel ==============================
B = zeros(Ntot,1);
trigIdx = zeros(C, nF);                                 % 0-based rising edges
for m = 1:nF
    for c = 1:C
        t0 = sm(m) + (c-1)*R + Dg;
        trigIdx(c,m) = t0;
        B(t0 + (1:Cfg.TrigWidth)) = Cfg.TrigAmp_V;
    end
end
if Cfg.TrigWidth >= R
    error('SBM:Trig', 'TrigWidth must be shorter than the record (%d).', R);
end

%% ======================= Excitation channel =============================
n = (0:Ntot-1)';
E = []; excIdx = [];
switch mode
    case "sine"
        nr  = Cfg.ExcRampPeriods*P;
        env = ones(Ntot,1);
        if nr > 0
            env(1:nr) = (1 - cos(pi*(0:nr-1)'/nr))/2;
            offStart = ceil(lastAcq/P)*P + P;           % keep ON through the last record
            k2 = offStart + (0:nr-1)';
            k2 = k2(k2 < Ntot);
            env(k2+1) = (1 + cos(pi*(0:numel(k2)-1)'/nr))/2;
            env(offStart+numel(k2)+1:end) = 0;
        end
        if s1 < nr + 5*P
            warning('SBM:Preamble', 'Preamble shorter than ramp + 5 periods.');
        end
        % offset faded with the sine: ao3 starts and ends at 0 V (no step on the piezo)
        E = (Cfg.ExcOffset_V + Cfg.ExcVpp_V/2 * sin(2*pi*n/P)) .* env;
    case "ttl"
        E = zeros(Ntot,1);
        excIdx = sm - Lead - (0:nF-1)'*Nphi;
        for m = 1:nF
            E(excIdx(m) + (1:Cfg.TTL_Width)) = Cfg.TTL_Amp_V;
        end
    case "gen_ao3"
        % external generator in 1-cycle triggered burst mode, re-triggered
        % every P A-lines by a TTL on ao3 -> sine exactly P-periodic on the
        % A-line clock (generator clock accuracy irrelevant)
        E = zeros(Ntot,1);
        excIdx = (0:P:(lastAcq + 2*P))';
        excIdx = excIdx(excIdx + Cfg.TTL_Width < Ntot);
        for k = 1:numel(excIdx)
            E(excIdx(k) + (1:Cfg.TTL_Width)) = Cfg.TTL_Amp_V;
        end
    case "gen_puff"
        % same generator setup; the ÷P TTL on ao3 is LabVIEW's repeated
        % micro-puff train (period P A-lines from sample PuffDelay_Alines,
        % over the whole waveform) -> no 4th row here, LabVIEW appends it
        d = 0;  if isfield(Cfg, 'PuffDelay_Alines'), d = Cfg.PuffDelay_Alines; end
        excIdx = (d:P:(lastAcq + 2*P))';
    case {"gen_counter", "none"}
        % gen_counter: same generator setup, but the ÷P trigger comes from a
        % DAQ counter clocked by the A-line trigger (no 4th AO row)
    otherwise
        error('SBM:ExcMode', ['ExcMode must be "gen_puff", "gen_counter", "gen_ao3", "sine", ' ...
            '"ttl" or "none".']);
end

waveform = [X'; Y'; B'];
if ~isempty(E), waveform = [waveform; E']; end
if any(abs(waveform(:)) > 10)
    error('SBM:Range', 'AO waveform exceeds +-10 V.');
end

%% ======================= Pixel / time bookkeeping =======================
% Position j (0-based) inside the frame body is acquired in record
% c = floor(j/R), A-line k = mod(j,R), i.e. at global index
%   n = s_m + j + Dg + lat   (beam lags the command by ~Dg, so the beam sits
%                              at body sample j when A-line n is recorded).
jAl = find(isAl) - 1;
PixelIdx = zeros(nB, nA);                               % 1-based within frame data
PixelIdx(sub2ind([nB nA], rowOf(isAl), colOf(isAl))) = jAl + 1;
tau = (PixelIdx - 1 + Dg + lat);                        % samples after s_m
Info = struct();
Info.Cfg = Cfg;
Info.P = P;  Info.Nphi = Nphi;  Info.S = S;  Info.Nbody = Nbody;
Info.G = G;  Info.R = R;  Info.C = C;  Info.K = K;  Info.F = F;
Info.FrameStart = sm;                                   % 0-based
Info.TrigIdx = trigIdx;  Info.ExcTrigIdx = excIdx;
Info.BurstCycles = BurstCycles;  Info.LeadSamples = Lead;
Info.ExcMode = mode;
Info.ExcPeriodic = ismember(mode, ["sine", "gen_ao3", "gen_counter", "gen_puff"]);  % P-periodic on the A-line clock
Info.Ntot = Ntot;  Info.Duration_s = Ntot*Ts;
Info.AcqWindow_s = (lastAcq - s1)*Ts;
Info.f_exc = 1/Cfg.T_excitation;
Info.Ts_SBM = Nphi*Ts;
Info.PixelIdx = PixelIdx;                               % (b, col) -> sample in frame record data
Info.TauSamples = tau;                                  % acquisition delay of each pixel in frame
ph0 = 0;                                                % excitation phase origin (trigger sample)
if mode == "gen_puff" && isfield(Cfg, 'PuffDelay_Alines'), ph0 = Cfg.PuffDelay_Alines; end
Info.PhiAcq = mod(2*pi*(s1 + tau - ph0)/P + pi, 2*pi) - pi;   % accumulated acq. phase (global ref.)
Info.FramePhaseSamples = mod(s1 + (0:nF-1)'*Nphi, P);
Info.xAxis_V = linspace(Vx(1), Vx(2), nA);
Info.yAxis_V = linspace(Vy(2), Vy(1), nB);              % row 1 = top
Info.IsAlineBody = isAl;
Info.RunInSamples = Nd;  Info.FlybackSamples = Lfly;
Info.Overshoot_V = overshoot .* [diff(Vx) diff(Vy)]/2;  % beyond the FOV edge [V]
Info.TotalAlines = R*C*nF;                              % -> TotCluster_Scans
Info.LabVIEW = struct( ...
    'A_lines_per_B_frame', R, ...                       % record length per B-trigger
    'B_frames', C, ...                                  % records per SBM frame
    'No_Cluster', nF, ...                               % SBM frames
    'No_Sync_Samples', nT, ...
    'B_frame_Signal_Delay', Dg, ...
    'TotCluster_Scans', R*C*nF, ...
    'AO_channels', size(waveform,1));
end

%% ========================================================================
function x = quintic(s, p0, v0, p1, v1)
% Quintic Hermite with zero end accelerations; v0, v1 in units of dx/ds.
s3 = s.^3;  s4 = s.^4;  s5 = s.^5;
x = (1 - 10*s3 + 15*s4 - 6*s5)*p0 + (s - 6*s3 + 8*s4 - 3*s5)*v0 ...
  + (-4*s3 + 7*s4 - 3*s5)*v1 + (10*s3 - 15*s4 + 6*s5)*p1;
end

function N = feasibleFrames(P, NP)
c = 1:(NP*P);
N = c(mod(NP*P, c) == 0 & (NP*P)./c < P/2);
end
