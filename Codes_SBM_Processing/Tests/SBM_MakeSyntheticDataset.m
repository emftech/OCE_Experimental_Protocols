function List = SBM_MakeSyntheticDataset(ExpDir, Exps, protocolDir, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_MakeSyntheticDataset
%
% Writes synthetic SBM acquisitions in the exact LabVIEW .bin format (text
% header + uint16 fringes, R x C x N_F A-lines), played through the real
% AO timeline of SBM_TwinOCT_Waveform (galvo lag = B-frame Signal Delay),
% to test the processing chain with a known ground truth:
%   - tilted surface + a weaker deeper layer, static speckle phase
%   - cylindrical wave from a needle at Source_mm, c = c0 (default 10.54 m/s,
%     the phantom's c_R), amplitude decay, + a reflection from a boundary
%   - fringe noise and a slow bulk drift
% Sidecars are stored three ways (to test SBM_FindSidecar):
%   'next'    -> Data/<name>_Sidecar.mat
%   'stamped' -> Sidecars/SBM_Waveform_<acquisition time - 20 s>.mat
%   'last'    -> Sidecars/SBM_Waveform_LAST.mat
%
%   Exps = struct('Period_us',{285,220}, 'NB',{21,21}, 'Rep',{1,1}, 'Sidecar',{'next','stamped'});
%   SBM_MakeSyntheticDataset('/tmp/SBM_Experiments', Exps, protocolDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4, opts = struct(); end
def = struct('nS', 512, 'c0', 10.54, 'Source_mm', [0.7 1.1], 'Amp_nm', 20, ...
    'Refl', 0.3, 'Wall_mm', 7, 'Noise', 25, 'Drift_nm', 10, 'FOV_mm', 12, 'NA', 37, ...
    't0', datetime(2026, 10, 18, 23, 30, 0), 'Seed', 3);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end, end
addpath(protocolDir);
dataDir = fullfile(ExpDir, 'Data');  scDir = fullfile(ExpDir, 'Sidecars');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
if ~exist(scDir, 'dir'), mkdir(scDir); end
rng(opts.Seed);
nS = opts.nS;  lam = 1.3e-6;  nidx = 1.39;
List = {};
warning('off', 'SBM:BscanPhase');
for e = 1:numel(Exps)
    X = Exps(e);
    Cfg = SBM_TwinOCT_Config();
    T = X.Period_us;  P = T/5;
    k = 1; while k*P - opts.NA < 20, k = k + 1; end
    Cfg.T_excitation = T*1e-6;  Cfg.nAscans = opts.NA;  Cfg.nTurnSamples = k*P - opts.NA;
    Cfg.nBscans = X.NB;  Cfg.nFrames = P;  Cfg.nPeriods = 4;
    Cfg.horizontalMin_V = -opts.FOV_mm/2;  Cfg.horizontalMax_V = opts.FOV_mm/2;
    Cfg.verticalMin_V = -opts.FOV_mm/2;  Cfg.verticalMax_V = opts.FOV_mm/2;
    Cfg.HVC = 1;  Cfg.VVC = 1;
    [W, Info] = SBM_TwinOCT_Waveform(Cfg);
    waveform = W;
    name = sprintf('PC_H_%dus_1c_1Vpp_%d_%d_%d_SBM_%d_%d_%d_%d_r%d', T, Cfg.nAscans, ...
        Cfg.nTurnSamples, Cfg.nBscans, Cfg.nFrames, Cfg.TrigDelay, opts.FOV_mm, opts.FOV_mm, X.Rep);
    tAcq = opts.t0 + minutes(2*(e-1));
    R = Info.R;  C = Info.C;  nF = Cfg.nFrames;  Ts = Cfg.Ts_OCT;
    w0 = 2*pi/(T*1e-6);  k0 = w0/(opts.c0*1e3);         % rad/mm

    keys = {'Header length in bytes (1 char <-> 1 byte)', 'Date/Time scan acquisition start', ...
        'Sample trigger source', 'Sample Rate (MHz)', '# spectral samples per A-line', ...
        '# pre-trigger samples', 'Scan pattern control', 'Scan type', 'A-lines per B-frame', ...
        'Flyback duration in # A-lines', 'B-frame trigger delay in # A-lines', ...
        '# of B-frames in 3D scan', '# of 3D repeated scans', 'Horizontal scan length (mm)', ...
        'Vertical scan length (mm)', 'Max A-lines per buffer'};
    vals = {'0', char(datetime(tAcq, 'Format', 'M/d/yyyy h:mm:ss a', 'Locale', 'en_US')), ...
        'Fast External', '1000.000', num2str(nS), '0', 'GUI-generated scan pattern waveform', ...
        'Linear', num2str(R), num2str(Cfg.nTurnSamples), num2str(Cfg.TrigDelay), num2str(C), ...
        num2str(nF), sprintf('%.3f', opts.FOV_mm), sprintf('%.3f', opts.FOV_mm), '1446'};
    H = '';
    for i = 1:numel(keys), H = [H keys{i} ': ' vals{i} ';' char(13)]; end %#ok<AGROW>
    H = strrep(H, 'Header length in bytes', [sprintf('%d', numel(H)+5) ': Header length in bytes']);
    fid = fopen(fullfile(dataDir, [name '.bin']), 'w');
    fwrite(fid, numel(H), 'uint16');  fwrite(fid, 0, 'uint16');  fwrite(fid, H, 'char');
    kk = (0:nS-1)';
    spk1 = 2*pi*rand(1, R*C);  spk2 = 2*pi*rand(1, R*C);
    for mfr = 1:nF
        for c = 1:C
            n = Info.TrigIdx(c, mfr) + Cfg.AlazarLatency + (0:R-1);
            cmd = max(n - Cfg.TrigDelay - Cfg.AlazarLatency, 0);
            x = W(1, cmd+1);  y = W(2, cmd+1);
            r  = hypot(x - opts.Source_mm(1), y - opts.Source_mm(2));
            rm = hypot(x - (2*opts.Wall_mm - opts.Source_mm(1)), y - opts.Source_mm(2));
            u = opts.Amp_nm*1e-9*(exp(-r/8).*sin(w0*n*Ts - k0*r) + ...
                opts.Refl*exp(-rm/8).*sin(w0*n*Ts - k0*rm)) + ...
                opts.Drift_nm*1e-9*sin(2*pi*1.3*n*Ts);
            ph = 4*pi*nidx*u/lam;
            z0 = 110 + 10*x/6 + 4*y/6;
            refl = 600*(1 + 0.3*cos(2*pi*x/4));
            sel = (c-1)*R + (1:R);
            fr = 2000 + refl.*cos(2*pi*kk*z0/nS + ph + spk1(sel)) + ...
                 0.4*refl.*cos(2*pi*kk*(z0+40)/nS + ph + spk2(sel)) + opts.Noise*randn(nS, R);
            fwrite(fid, uint16(round(fr)), 'uint16');
        end
    end
    fclose(fid);
    switch X.Sidecar
        case 'next'
            save(fullfile(dataDir, [name '_Sidecar.mat']), 'Cfg', 'Info', 'waveform', '-v7.3');
        case 'stamped'
            save(fullfile(scDir, ['SBM_Waveform_' char(datetime(tAcq - seconds(20), 'Format', 'yyyyMMdd_HHmmss')) '.mat']), ...
                'Cfg', 'Info', 'waveform', '-v7.3');
        case 'last'
            save(fullfile(scDir, 'SBM_Waveform_LAST.mat'), 'Cfg', 'Info', 'waveform', '-v7.3');
    end
    List{end+1} = name; %#ok<AGROW>
    fprintf('wrote %s (%s sidecar)\n', name, X.Sidecar);
end
warning('on', 'SBM:BscanPhase');
end
