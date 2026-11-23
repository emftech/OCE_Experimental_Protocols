function G = SBM_TwinOCT_GeneratorSettings(Cfg, fAline_Hz, model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 18-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_GeneratorSettings
%
% Prints the exact settings for the external function generator and for
% the ÷P counter trigger, for the excitation modes that use a generator:
%
%   "gen_counter" / "gen_ao3"  (recommended): the generator is in TRIGGERED
%       BURST mode with N = 1 cycle, re-triggered every P A-lines. Each cycle
%       restarts on the A-line clock, so the excitation is exactly P-periodic
%       whatever the generator's own clock accuracy. The generator is set by
%       PERIOD, slightly SHORTER than the trigger period T_trig = P/f_Aline:
%           T_gen = T_trig - Cfg.Gen_DeadTime_s   (282.000 us at 200 kHz)
%       so every cycle ends (and the generator re-arms) before the next
%       trigger; the output holds the start level (0 V) for the dead time.
%   "ttl": one trigger per frame, N-cycle burst at f_trig (free-running
%       inside a frame; frequency must match the A-line clock closely).
%
%   G = SBM_TwinOCT_GeneratorSettings                 % config, 200 kHz nominal
%   G = SBM_TwinOCT_GeneratorSettings(Cfg, 199912.4, "DG922Pro")  % + Rigol DG900 Pro
%       key sequence (the lab's generator)
%   G = SBM_TwinOCT_GeneratorSettings(Cfg, 199912.4)  % with the MEASURED A-line rate
%       (Out.AlineRateFromTrig_Hz of SBM_TwinOCT_CheckExcitation, or the laser's
%        sweep-rate readout)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1 || isempty(Cfg), Cfg = SBM_TwinOCT_Config(); end
if nargin < 2 || isempty(fAline_Hz), fAline_Hz = 1/Cfg.Ts_OCT; end
if nargin < 3, model = ""; end
warning('off', 'SBM:BscanPhase');
[~, Info] = SBM_TwinOCT_Waveform(Cfg);
warning('on', 'SBM:BscanPhase');
P = Info.P;  mode = lower(string(Cfg.ExcMode));
fTrig = fAline_Hz / P;  Ttrig = 1/fTrig;

G = struct();
G.Mode = mode;  G.AlineRate_Hz = fAline_Hz;  G.P = P;
G.TriggerRate_Hz = fTrig;
G.Amplitude_Vpp = Cfg.ExcVpp_V;  G.Offset_V = Cfg.ExcOffset_V;

fprintf('\n=========== SBM external generator settings (%s) ===========\n', mode);
fprintf('A-line clock %.3f Hz, P = %d  ->  excitation period %.4f us (f_exc = %.4f Hz)\n', ...
    fAline_Hz, P, 1e6*Ttrig, fTrig);

switch mode
    case {"gen_counter", "gen_ao3", "gen_puff"}
        dead = deadTime(Cfg);
        Tgen = Ttrig - dead;
        fGen = 1/Tgen;
        G.Period_us = 1e6*Tgen;
        G.Frequency_Hz = fGen;
        G.Cycles = 1;
        G.DeadTime_us = 1e6*dead;
        G.RearmOK = dead > Cfg.Gen_RearmTime_s;
        % distortion caused by the dead time (one exact period, fine grid)
        Nfine = 2^14;  tt = (0:Nfine-1)'/Nfine*Ttrig;
        x = sin(2*pi*fGen*tt) .* (tt < 1/fGen);
        Xf = fft(x)/Nfine;
        a1 = 2*abs(Xf(2));
        G.FundamentalVpp_per_Vpp = a1;                   % amplitude of the 1st harmonic
        G.THD_dB = 20*log10(sqrt(sum(abs(2*Xf(3:Nfine/2)).^2))/a1);
        fprintf('\n--- Function generator ---\n');
        fprintf('  Function            : SINE\n');
        fprintf('  PERIOD (set this)   : %.3f us   (= %.3f us trigger period - %.3f us dead time)\n', ...
            1e6*Tgen, 1e6*Ttrig, 1e6*dead);
        fprintf('  (equivalent freq.   : %.4f Hz, i.e. f_exc x %.4f)\n', fGen, fGen/fTrig);
        fprintf('  Amplitude           : %.3f Vpp   (fundamental delivered: %.3f Vpp)\n', ...
            Cfg.ExcVpp_V, Cfg.ExcVpp_V*a1);
        fprintf('  Offset              : %.3f V\n', Cfg.ExcOffset_V);
        fprintf('  Output load         : HIGH-Z (amplifier input), not 50 ohm\n');
        fprintf('  Burst               : ON, mode TRIGGERED (N-cycle), cycles = 1\n');
        fprintf('  Burst start phase   : 0 deg\n');
        fprintf('  Trigger source      : EXTERNAL, slope RISING, trigger delay 0\n');
        fprintf('  Idle level          : start phase / 0 V (if the menu offers it)\n');
        fprintf('  Dead time per period: %.2f us  (re-arm needs %.2f us)  [%s]\n', 1e6*dead, ...
            1e6*Cfg.Gen_RearmTime_s, okStr(G.RearmOK));
        fprintf('  Distortion from dead time: THD %.1f dB\n', G.THD_dB);
        if ~G.RearmOK
            fprintf(2, '  -> increase Cfg.Gen_DeadTime_s (e.g. 5e-6) or the generator skips cycles.\n');
        end
        G.CounterHighTicks = Cfg.TTL_Width;  G.CounterLowTicks = P - Cfg.TTL_Width;
        if mode == "gen_counter"
            fprintf('\n--- DAQ counter trigger (SBM_CounterTrigger VI) ---\n');
            fprintf('  Channel             : Dev1/ctr0 (any free counter)  CO Pulse Ticks\n');
            fprintf('  Source of ticks     : the PFI terminal that carries the A-line clock\n');
            fprintf('                        (the same terminal used as the AO sample clock)\n');
            fprintf('  High ticks / Low ticks: %d / %d   (period = %d A-lines)\n', ...
                G.CounterHighTicks, G.CounterLowTicks, P);
            fprintf('  Idle state          : LOW, initial delay 0 ticks\n');
            fprintf('  Timing              : Implicit, CONTINUOUS samples\n');
            fprintf('  Output              : the counter''s output terminal (e.g. Dev1/Ctr0Out) -> generator EXT TRIG\n');
        elseif mode == "gen_puff"
            d = 0;  if isfield(Cfg, 'PuffDelay_Alines'), d = Cfg.PuffDelay_Alines; end
            G.PuffFrequency_Hz = fTrig;
            G.PuffOnLength_ms = 1e3*Cfg.TTL_Width/fAline_Hz;
            G.PuffDelay_ms = 1e3*d/fAline_Hz;
            fprintf('\n--- LabVIEW Air puff controller (master front panel) ---\n');
            fprintf('  Enable Puff         : ON\n');
            fprintf('  Enum (Transducer)   : Repeated micro puffs\n');
            fprintf('  Puff Frequency (Hz) : %.4f   (= f_Aline/%d; repetition period must come out = %d A-lines)\n', fTrig, P, P);
            fprintf('  Puff ON Length (ms) : %.3f    (%d A-lines, < period)\n', G.PuffOnLength_ms, Cfg.TTL_Width);
            fprintf('  Puff Delay (ms)     : %.3f    (%d A-lines)\n', G.PuffDelay_ms, d);
            fprintf('  Sweep start/end/step, bin length: 0 (unused). Macro puff voltage: as in QH (TTL level).\n');
            fprintf('  Scope check: ao3 = continuous %.1f us pulse train for the WHOLE scan (~%.0f ms),\n', 1e6*Ttrig, 1e3*Info.Duration_s);
            fprintf('               generator output = one sine cycle per pulse (no gaps, no 290/280 us periods).\n');
        else
            fprintf('\n--- ao3 trigger (row 4 of the SBM waveform) ---\n');
            fprintf('  %.0f V pulses, %d A-lines wide, every %d A-lines -> generator EXT TRIG\n', ...
                Cfg.TTL_Amp_V, Cfg.TTL_Width, P);
        end
    case "ttl"
        G.Frequency_Hz = fTrig;  G.Cycles = Info.BurstCycles;
        fprintf('\n--- Function generator (per-frame burst, fallback) ---\n');
        fprintf('  SINE, %.4f Hz, %.3f Vpp, HIGH-Z, Burst TRIGGERED, cycles = %d, phase 0, EXT rising\n', ...
            fTrig, Cfg.ExcVpp_V, Info.BurstCycles);
        fprintf('  Frequency must match f_Aline/P: a 100 ppm error drifts %.2f deg over one frame.\n', ...
            360*1e-4*Info.BurstCycles);
    otherwise
        fprintf('ExcMode "%s" does not use an external generator.\n', mode);
end
if contains(lower(string(model)), ["dg9", "dg8"]) && isfield(G, 'Cycles') && G.Cycles == 1
    fprintf('\n--- Rigol DG900 Pro / DG800 Pro sequence (CH1) ---\n');
    fprintf('  1. [Default] -> confirm. CH1 tile: Imped = HighZ (keep).\n');
    fprintf('  2. [Sine]. Tap the Freq/Period key to show PERIOD, tap it -> keypad %.3f -> unit us.\n', 1e6*Tgen);
    fprintf('  3. Tap Ampl -> %.3f -> Vpp.  Offset -> %.3f Vdc.  Phase -> 0 deg.\n', Cfg.ExcVpp_V, Cfg.ExcOffset_V);
    fprintf('  4. Mode dropdown "Cont" -> Burst  (or press [Burst]).\n');
    fprintf('  5. Burst: Type = N Cycle, Cycles = 1, Burst/Start Phase = 0 deg,\n');
    fprintf('     Trigger Source = External, Edge = Rising (Rise), Delay = 0 s,\n');
    fprintf('     Idle Level = 1st Pt (the output rests at the start of the cycle between bursts).\n');
    fprintf('  6. ÷%d TTL (ao3) -> the CH1 EXTERNAL TRIGGER input used in QH MB mode (same cable);\n', P);
    fprintf('     NOT the front Counter BNC. CH1 output -> amplifier input (as in QH).\n');
    fprintf('  7. [1 On/Off] -> output ON.  [Store/Recall] -> save as "SBM_%dus".\n', round(1e6*Cfg.T_excitation));
    if mode == "gen_ao3"
        fprintf('  Check: front Counter input (1 MOhm, DC) on the A-line trigger (PFI0) reads f_Aline;\n');
        fprintf('         re-run with that value to get the exact period.\n');
    else
        fprintf('  Check: front Counter input (1 MOhm, DC) on the ÷%d trigger must read %.4f Hz.\n', P, fTrig);
    end
end
fprintf('=============================================================\n');
end

function d = deadTime(Cfg)
% dead time per period; older configs used a frequency margin instead
if isfield(Cfg, 'Gen_DeadTime_s') && ~isempty(Cfg.Gen_DeadTime_s)
    d = Cfg.Gen_DeadTime_s;
elseif isfield(Cfg, 'Gen_FreqMargin')
    d = Cfg.T_excitation*(1 - 1/(1 + Cfg.Gen_FreqMargin));
else
    d = 3e-6;
end
end

function s = okStr(b)
if b, s = 'OK'; else, s = 'TOO SHORT'; end
end
