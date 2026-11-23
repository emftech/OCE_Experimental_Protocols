%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_SmokeTestNull
%
% Tests SBM_TwinOCT_NullCheck on a synthetic NO-EXCITATION acquisition in
% the exact LabVIEW .bin format: tilted surface with two bumps and a
% striped reflectivity, static speckle phase, finite SNR, 20 nm bulk drift,
% and a galvo that lags the command by B-frame Signal Delay + 3 A-lines.
% Expected: the check recommends delay = 33 (configured 30), and reports
% a displacement noise consistent with the SNR.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath'));  addpath(fileparts(here));  % SBM_Protocol root
out = fullfile(tempdir, 'sbm_null_smoke');  if ~exist(out,'dir'), mkdir(out); end
Cfg = SBM_TwinOCT_Config();
Cfg.horizontalMin_V = -6; Cfg.horizontalMax_V = 6; Cfg.verticalMin_V = -6; Cfg.verticalMax_V = 6;
Cfg.ExcMode = "none";
[W, Info] = SBM_TwinOCT_Waveform(Cfg);
waveform = W;  save(fullfile(out, 'null_Sidecar.mat'), 'Cfg', 'Info', 'waveform');
nS = 512;  R = Info.R;  C = Info.C;  nF = Cfg.nFrames;  Ts = Cfg.Ts_OCT;
lagTrue = Cfg.TrigDelay + 3;                           % galvo slower than assumed
lam = 1.3e-6; nidx = 1.39;

keys = {'Header length in bytes (1 char <-> 1 byte)','Date/Time scan acquisition start', ...
    'Sample trigger source','Sample Rate (MHz)','# spectral samples per A-line', ...
    '# pre-trigger samples','Scan pattern control','Scan type','A-lines per B-frame', ...
    'Flyback duration in # A-lines','B-frame trigger delay in # A-lines', ...
    '# of B-frames in 3D scan','# of 3D repeated scans','Horizontal scan length (mm)', ...
    'Vertical scan length (mm)'};
vals = {'0','10/17/2026 6:30:00 PM','Fast External','1000.000',num2str(nS),'0', ...
    'GUI-generated scan pattern waveform','Linear',num2str(R),num2str(Cfg.nTurnSamples), ...
    num2str(Cfg.TrigDelay),num2str(C),num2str(nF),'12.000','12.000'};
H = '';  for i = 1:numel(keys), H = [H keys{i} ': ' vals{i} ';' char(13)]; end %#ok<AGROW>
H = strrep(H, 'Header length in bytes', [sprintf('%d', numel(H)+5) ': Header length in bytes']);
fid = fopen(fullfile(out, 'null.bin'), 'w');
fwrite(fid, numel(H), 'uint16');  fwrite(fid, 0, 'uint16');  fwrite(fid, H, 'char');

rng(2);
kk = (0:nS-1)';
zSurf = @(x,y) 110 + 15*x/6 + 10*exp(-((x-2).^2+(y+1).^2)/2) - 8*exp(-((x+3).^2+(y-3).^2)/1.5);
refl  = @(x,y) 500*(1 + 0.6*cos(2*pi*x/3)) .* (1 + 0.3*(y > 0));
specklePh = 2*pi*rand(1, R*C);
for m = 1:nF
    for c = 1:C
        n   = Info.TrigIdx(c,m) + Cfg.AlazarLatency + (0:R-1);
        cmd = max(n - lagTrue - Cfg.AlazarLatency, 0);
        x = W(1, cmd+1);  y = W(2, cmd+1);
        z0 = zSurf(x, y);
        bulk = 4*pi*nidx*(20e-9*sin(2*pi*1.5*n*Ts))/lam;   % 20 nm, 1.5 Hz drift
        phi = specklePh((c-1)*R+(1:R)) + bulk;
        fr = 2000 + refl(x,y) .* cos(2*pi*kk*z0/nS + phi) + 25*randn(nS, R);
        fwrite(fid, uint16(round(fr)), 'uint16');
    end
end
fclose(fid);

Res = SBM_TwinOCT_NullCheck(fullfile(out, 'null.bin'), '', struct('Cut_Depth', [30 250]));
fprintf('\nSMOKE: recommended delay %d (expected %d)\n', Res.Summary.RecommendedDelay, lagTrue);
