function Cfg = SBM_TwinOCT_Config()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_Config
%
% Single source of truth for the Synchronized BM-mode (SBM) protocol on the
% Twin OCT (1300 nm SS-OCT, ATS9360 + NI-DAQ AO clocked by the A-line
% trigger). Edit THIS file to change the protocol. Both the LabVIEW MATLAB
% script node (Galvo 3D cluster driver_SBM.vi) and the offline tools read it.
%
% Thesis notation (Chapter 5):
%   P      = T_excitation / Ts_OCT   samples per excitation period (57)
%   N_F    = nFrames                 synchronized BM frames (57)
%   N_P    = nPeriods                reconstructed excitation periods (4)
%   N_phi  = N_P * P / N_F           phase-step samples (4)  -> Eq. (46)
%   N_A    = nAscans, N_B = nBscans, N_turn = nTurnSamples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------------------------- Clock & excitation ------------------------
Cfg.Ts_OCT       = 5e-6;      % A-line period [s] (200 kHz sweep = DAQ AO clock)
Cfg.T_excitation = 285e-6;    % excitation period [s] -> P = 57 samples

%% ---------------------------- SBM sampling ------------------------------
Cfg.nAscans      = 37;        % valid A-lines per B-scan (N_A)
Cfg.nTurnSamples = 20;        % turn samples between B-scans (N_turn)
                              % N_A + N_turn = 57 = P -> every B-scan starts
                              % at the same excitation phase (thesis Sec. 5.4)
Cfg.nBscans      = 21;        % B-scans per BM frame (N_B): 11/21/41/81
Cfg.nFrames      = 57;        % synchronized BM frames (N_F)
Cfg.nPeriods     = 4;         % reconstructed periods (N_P) -> N_phi = 4

%% ---------------------------- Field of view -----------------------------
% Volts at the galvo inputs. In LabVIEW these are overwritten by the
% "Horizontal/Vertical Min/Max (V)" values that XY_Scanning_Rectangle_SDK.vi
% computes from the scan length (mm) and HVC/VVC (V/mm).
Cfg.horizontalMin_V = -1;
Cfg.horizontalMax_V =  1;
Cfg.verticalMin_V   = -1;
Cfg.verticalMax_V   =  1;
Cfg.HVC = [];                 % optional V/mm (only for mm axes in processing)
Cfg.VVC = [];

%% ---------------------------- Alazar records ----------------------------
% The ATS9360 runs in AUX_IN_TRIGGER_ENABLE mode: every rising edge on the
% B-trigger line (DAQ AO "ref" row) arms one record of "A-lines per
% B-frame" consecutive A-lines. SBM groups G consecutive B-scan steps per
% record: R = G*(N_A+N_turn) A-lines per trigger, C = ceil(N_B/G) triggers
% per frame. One record must fit in one DMA buffer (8 MB):
Cfg.MaxAlinesPerBuffer = 1446; % from the 06-10-2026 headers (2880 samples/A-line)
Cfg.BscansPerRecord    = [];   % [] = automatic (as many as fit in a buffer)

%% ---------------------------- Triggers ----------------------------------
Cfg.TrigDelay   = 30;         % "B-frame Signal Delay" [A-lines]: galvo lag
Cfg.TrigWidth   = 10;         % B-trigger pulse width [A-lines]
Cfg.TrigAmp_V   = 5;          % B-trigger level [V] (same as MB "ref" = 5 V)
Cfg.AlazarLatency = 0;        % A-lines between AUX edge and first record
                              % (constant -> only a global phase; measure it)

%% ---------------------------- Excitation ------------------------------
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>  CHOOSE THE VARIANT HERE  <<<<<<<<<<<<<<<<<<<<<<
% Two complete implementations live in this folder (22-11-2026):
%   "external_generator" (Variant_A_ExternalGenerator, DEFAULT now):
%       Rigol DG922 Pro, sine set by PERIOD (P/f_Aline - 3 us = 282 us),
%       1-cycle N-Cycle burst, EXTERNAL rising trigger every P A-lines.
%       The trigger source is chosen with Cfg.GenTrigger:
%         "puff"    -> LabVIEW's own "Repeated micro puffs" pulse train on ao3
%                      (Enable Puff ON, Puff Frequency = f_Aline/P). No LabVIEW
%                      block-diagram edit: same path as the QH MB mode. [default]
%         "ao3_row" -> TTL written as AO row 4 of the SBM waveform; needs the
%                      4-channel AO task (case "AsyncOCE","SBM" in
%                      LabVIEW_Setup_Mirrors_and_Btrig.vi, Enable Puff OFF).
%         "counter" -> ÷P TTL from a DAQ counter (SBM_CounterTrigger.vi).
%   "daq_sine" (Variant_B_DAQ_Sine): the NI DAQ synthesizes the sine itself on
%       ao3 (AsyncOCE-style, 4 AO rows), no generator; needs the same
%       "AsyncOCE","SBM" case edit. Steps: Variant_B_DAQ_Sine/*.md.
% The processing (Codes_SBM_Processing) is the same for both variants.
Cfg.Variant    = "external_generator";   % "external_generator" | "daq_sine"
Cfg.GenTrigger = "puff";                 % "puff" | "ao3_row" | "counter" (external generator only)
%
% Low-level modes (set automatically from Cfg.Variant at the end of this
% file; a script can still override Cfg.ExcMode afterwards, e.g. "none"):
% SBM needs the excitation to be EXACTLY periodic with P = 57 A-lines.
% "sine"   : (24-10-2026, AsyncOCE-style). The DAQ synthesizes the
%         sine itself on AO row 4 = ao3, sampled by the A-line clock (57
%         samples per period at 285 us), so it is P-periodic by construction
%         (no external generator, no air-puff controller). ao3 -> amplifier
%         input (the cable that went to the Rigol OUTPUT side). Needs the
%         4-channel AO task in LabVIEW_Setup_Mirrors_and_Btrig.vi (the
%         "AsyncOCE" case, edited to also accept "SBM"; Enable Puff OFF).
%         Level: ExcOffset_V + ExcVpp_V/2 sin(2 pi n/P), the whole signal
%         (offset included) faded in/out over ExcRampPeriods, so ao3 starts
%         and ends at 0 V. Default 0-1 V = the Rigol drive of the QH runs
%         (1 Vpp, 500 mV offset); the old AsyncOCE driver used 0-0.95 V.
% "gen_puff": (20-10-2026) Same generator setup (1-cycle triggered
%         burst, PERIOD = P/f_Aline - dead time), but the ÷P TTL on ao3 is
%         built by LabVIEW's OWN air-puff logic (LabVIEW_Setup_Mirrors_and_Btrig.vi,
%         "Repeated micro puffs", Enable Puff ON): a valve-ON pulse every
%         "puff repetition period in multiples of A-lines" over the whole
%         galvo waveform, on the same AO task/clock as the galvos. This is the
%         path the QH MB mode already uses (proven on the scope 20-10-2026).
%         The SBM waveform has 3 rows; LabVIEW appends the puff row as ao3.
%         Panel: Puff Frequency = f_Aline/P (3508.7719 Hz), Puff ON Length =
%         TTL_Width/f_Aline (0.05 ms), Puff Delay = PuffDelay_Alines/f_Aline.
%         Print them with SBM_TwinOCT_GeneratorSettings.
% "gen_ao3": (19-10-2026; ao3 stayed silent on the Twin) External function generator in
%         triggered 1-cycle burst mode, re-triggered every P A-lines by a
%         5 V TTL on AO row 4 = ao3 (the same BNC that triggers the generator
%         in the QH multimeridian MB mode). Needs the 4-channel AO task: in
%         the Twin's LabVIEW that is Air-puff controller Transducer =
%         "AsyncOCE" (Enable Puff OFF). The trigger only exists while the AO
%         buffer plays, so the excitation starts PreamblePeriods before frame 1
%         and stops by itself after the acquisition.
% "gen_counter": same generator setup, ÷P TTL from a DAQ counter clocked by
%         the A-line trigger (SBM_CounterTrigger.vi); 3-row waveform.
% "sine"   : the DAQ synthesizes the sine on ao3 (AsyncOCE-style).
% "ttl"    : one TTL per frame, generator in N-cycle burst (free-running
%            inside each frame; fallback only).
% "none"   : no excitation (null runs).
% NOTE: the k-clock CANNOT be used: it is the swept-source optical-frequency
% sample clock (hundreds of MHz, non-uniform in time, gated between
% sweeps). Lock to the A-line (sweep) trigger instead.
Cfg.ExcMode          = "";      % derived from Cfg.Variant / Cfg.GenTrigger below
Cfg.PuffDelay_Alines = 0;     % "gen_puff": LabVIEW Puff Delay [A-lines] (0 ms); only a
                              % global phase, but keep sidecar and panel equal
Cfg.ExcVpp_V         = 1.0;   % amplitude at the amplifier input [Vpp] (Rigol Ampl. / DAQ ao3)
Cfg.ExcOffset_V      = 0.5;   % offset [V]: 0.5 -> 0..1 V, the QH Rigol drive (Ampl 1 Vpp,
                              % Offset 500 mV); set 0 for a bipolar amplifier input
Cfg.Gen_DeadTime_s   = 3e-6;  % the generator is set by PERIOD:
                              %   T_gen = T_trig - Gen_DeadTime_s = 57/f_Aline - 3 us
                              %   (282.000 us at 200 kHz) so each 1-cycle burst ends
                              %   before the next trigger (re-arm time)
Cfg.Gen_RearmTime_s  = 1e-6;  % generator re-arm time after a burst (datasheet;
                              % 1 us is safe for most; check yours)
Cfg.ExcRampPeriods   = 5;     % "sine": Hann fade-in/out (avoids piezo step)
Cfg.PreamblePeriods  = 50;    % excitation ON before frame 1 (steady state; ~14 ms)
Cfg.PostPeriods      = 5;     % after the last record
Cfg.TTL_LeadPeriods  = 4;     % "ttl": periods between trigger and frame 1
Cfg.TTL_GapPeriods   = 2;     % "ttl": dead time between bursts
Cfg.TTL_Width        = 10;    % trigger pulse width [A-lines] ("ttl", "gen_ao3", counter)
Cfg.TTL_Amp_V        = 5;     % trigger pulse level [V]

%% ---------------------------- Galvo transitions -------------------------
% All galvo transitions are quintic (C2): position, velocity and
% acceleration are continuous at every junction (no kinks for the servo).
Cfg.FlybackMinSamples = 2*57; % rest-to-rest return between frames (min.)
Cfg.TurnReversalSamples = 21; % turns longer than this use stop-dwell-go (bounded overshoot)
Cfg.RunInSamples      = [];   % run-in/run-out at frame ends; [] = (min(N_turn,20)+1)/sqrt(2)
                              % (same peak jerk as a B-scan turn, lower acceleration)

%% ---------------------------- Bookkeeping -------------------------------
Cfg.SidecarDir = "";          % "" = no sidecar; LabVIEW node sets a folder

%% ---------------------------- Variant -> ExcMode ------------------------
switch lower(string(Cfg.Variant))
    case "external_generator"
        switch lower(string(Cfg.GenTrigger))
            case "puff",    Cfg.ExcMode = "gen_puff";
            case "ao3_row", Cfg.ExcMode = "gen_ao3";
            case "counter", Cfg.ExcMode = "gen_counter";
            otherwise, error('SBM:Variant', 'Cfg.GenTrigger must be "puff", "ao3_row" or "counter".');
        end
    case "daq_sine"
        Cfg.ExcMode = "sine";
    otherwise
        error('SBM:Variant', 'Cfg.Variant must be "external_generator" or "daq_sine".');
end

% make the tools / variant folders reachable once this folder is on the path
here = fileparts(mfilename('fullpath'));
for d = ["Tools" "Variant_A_ExternalGenerator" "Variant_B_DAQ_Sine"]
    if isfolder(fullfile(here, d)) && ~contains(path, fullfile(here, d))
        addpath(fullfile(here, d));
    end
end
end
