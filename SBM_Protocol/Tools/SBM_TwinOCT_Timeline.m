function T = SBM_TwinOCT_Timeline(Period_us, NA, NB, NP, NF, FOV_mm, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 20-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_Timeline
%
% Exact per-A-line sampling timeline of an SBM acquisition on the Twin OCT,
% for simulations. It is built by the SAME generator that drives the
% hardware (SBM_TwinOCT_Config + SBM_TwinOCT_Waveform), not by an idealized
% raster: quintic galvo transitions, stop-dwell-go turns, galvo lag
% (B-frame Signal Delay), records triggered by the B-trigger, frame period
% F = K P + N_phi. The beam position of A-line n is the galvo COMMAND at
% n - TrigDelay - AlazarLatency (the lag model used for phi_acq), as in
% Tests/SBM_MakeSyntheticDataset of Codes_SBM_Processing.
%
%   T = SBM_TwinOCT_Timeline(285, 37, 21, 4, 57, 12)
%   T = SBM_TwinOCT_Timeline(285, 37, 21, 4, 57, [18 18], struct('nTurn', 20))
%
% Inputs: excitation period [us], N_A, N_B, N_P, N_F, FOV [mm] (scalar or
%   [x y]); opts.nTurn (default: smallest k P - N_A >= 20, as the log),
%   opts.Cfg (start from this config instead of SBM_TwinOCT_Config()).
% Output struct T:
%   n        [Nrec x 1]  global sample index of every RECORDED A-line (0-based,
%                         sample 0 = first AO sample = excitation phase 0)
%   t_s      [Nrec x 1]  n*Ts [s]
%   x_mm, y_mm [Nrec x 1] beam position (scan centre = 0)
%   frame, rec, k        frame (1..N_F), record (1..C), A-line in record (0..R-1)
%   Pix.n    [N_B x N_A x N_F] sample index of every raster PIXEL, y ASCENDING
%            rows (same orientation as SBM_ToFrames / the processing)
%   Pix.x_mm [1 x N_A], Pix.y_mm [1 x N_B]  nominal pixel grid (ascending)
%   Info, Cfg, waveform : the sidecar exactly as the hardware would save it
%            (Info.PhiAcq = acquisition phase, Eq. 57; Info.PixelIdx, TauSamples)
% The excitation is P-periodic from sample 0: u(x,y,t) = A(r) cos(w n Ts - k r)
% sampled at T.Pix.n reproduces the real sampling; recover with
% SBM_HarmonicRecovery(Frames, Info, ...) of Codes_SBM_Processing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7, opts = struct(); end
if isscalar(FOV_mm), FOV_mm = [FOV_mm FOV_mm]; end
if isfield(opts, 'Cfg'), Cfg = opts.Cfg; else, Cfg = SBM_TwinOCT_Config(); end
Ts = Cfg.Ts_OCT;
P = round(Period_us*1e-6/Ts);
if isfield(opts, 'nTurn'), nT = opts.nTurn;
else, k = 1; while k*P - NA < 20, k = k + 1; end, nT = k*P - NA; end
Cfg.T_excitation = Period_us*1e-6;
Cfg.nAscans = NA;  Cfg.nTurnSamples = nT;  Cfg.nBscans = NB;
Cfg.nFrames = NF;  Cfg.nPeriods = NP;
% Galvos in "volts = mm x s": s keeps the waveform (with its ~42 % turn
% overshoot) inside the +-10 V check for any FOV; positions are rescaled.
s = 4 / max(FOV_mm/2);
Cfg.horizontalMin_V = -s*FOV_mm(1)/2;  Cfg.horizontalMax_V = s*FOV_mm(1)/2;
Cfg.verticalMin_V   = -s*FOV_mm(2)/2;  Cfg.verticalMax_V   = s*FOV_mm(2)/2;
Cfg.HVC = 1;  Cfg.VVC = 1;
ws = warning('off', 'SBM:BscanPhase');
[W, Info] = SBM_TwinOCT_Waveform(Cfg);
warning(ws);

R = Info.R;  C = Info.C;  lat = Cfg.AlazarLatency;  Dg = Cfg.TrigDelay;
[kk, cc, mm] = ndgrid(0:R-1, 1:C, 1:NF);
n = reshape(Info.TrigIdx(sub2ind([C NF], cc(:), mm(:))), [], 1) + lat + kk(:);
cmd = max(n - Dg - lat, 0) + 1;
cmd = min(cmd, size(W, 2));
T = struct();
T.n = n;  T.t_s = n*Ts;
T.x_mm = W(1, cmd)'/s;  T.y_mm = W(2, cmd)'/s;
T.frame = mm(:);  T.rec = cc(:);  T.k = kk(:);

% pixel instants (frame record index j = (c-1) R + k <-> Info.PixelIdx = j+1)
PIx = flipud(Info.PixelIdx);                    % y ascending
c0 = floor((PIx - 1)/R) + 1;  k0 = mod(PIx - 1, R);
Pn = zeros(NB, NA, NF);
for m = 1:NF
    Pn(:, :, m) = reshape(Info.TrigIdx(sub2ind([C NF], c0(:), m*ones(numel(c0), 1))), NB, NA) + lat + k0;
end
T.Pix.n = Pn;  T.Pix.t_s = Pn*Ts;
T.Pix.x_mm = linspace(-FOV_mm(1)/2, FOV_mm(1)/2, NA);
T.Pix.y_mm = linspace(-FOV_mm(2)/2, FOV_mm(2)/2, NB);
T.P = P;  T.Nphi = Info.Nphi;  T.F = Info.F;  T.Ts = Ts;  T.f_exc = 1/Cfg.T_excitation;
T.Info = Info;  T.Cfg = Cfg;  T.waveform = W;  T.VoltPerMm = s;
end
