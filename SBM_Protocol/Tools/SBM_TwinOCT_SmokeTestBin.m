%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_SmokeTestBin
%
% Writes a synthetic SBM acquisition in the exact LabVIEW .bin format
% (text header + uint16 fringes, R x C x N_F A-lines; one reflector with
% speckle-like static phase, displaced by an off-centre cylindrical wave
% played through the real AO timeline) into tempdir and runs
% SBM_TwinOCT_ProcessStarter on it. Expected: phase error ~0.01 rad and
% ~10 m/s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
here = fileparts(mfilename('fullpath'));  addpath(fileparts(here));  % SBM_Protocol root
out = fullfile(tempdir, 'sbm_smoke'); if ~exist(out,'dir'), mkdir(out); end
Cfg = SBM_TwinOCT_Config();
Cfg.horizontalMin_V=-6; Cfg.horizontalMax_V=6; Cfg.verticalMin_V=-6; Cfg.verticalMax_V=6; Cfg.HVC=1; Cfg.VVC=1;
[W, Info] = SBM_TwinOCT_Waveform(Cfg);
save(fullfile(out,'sidecar.mat'),'Cfg','Info','W');
nS = 512; R = Info.R; C = Info.C; nF = Cfg.nFrames; Ts = Cfg.Ts_OCT;
f0 = Info.f_exc; w0 = 2*pi*f0; c0 = 10; k0 = w0/(c0*1e3); src = [0.7 1.1];
lam = 1.3e-6; nidx = 1.39; uamp = 30e-9;
keys = {'Header length in bytes (1 char <-> 1 byte)','Date/Time scan acquisition start','Sample trigger source','Sample Rate (MHz)','# spectral samples per A-line','# pre-trigger samples','Scan pattern control','Scan type','A-lines per B-frame','Flyback duration in # A-lines','B-frame trigger delay in # A-lines','# of B-frames in 3D scan','# of 3D repeated scans','Horizontal scan length (mm)','Vertical scan length (mm)'};
vals = {'0','10/17/2026 10:00:00 AM','Fast External','1000.000',num2str(nS),'0','GUI-generated scan pattern waveform','SBM',num2str(R),num2str(Cfg.nTurnSamples),num2str(Cfg.TrigDelay),num2str(C),num2str(nF),'12.000','12.000'};
H = ''; for i=1:numel(keys), H = [H keys{i} ': ' vals{i} ';' char(13)]; end
H = strrep(H, 'Header length in bytes', [sprintf('%d',numel(H)+5) ': Header length in bytes']);
fid = fopen(fullfile(out,'smoke.bin'),'w'); fwrite(fid, numel(H), 'uint16'); fwrite(fid, 0, 'uint16'); fwrite(fid, H, 'char');
kk = (0:nS-1)'; z0 = 100; rng(1); phStatic = 2*pi*rand(1,R*C); % speckle-like static phase per raster position
for m = 1:nF
  for c = 1:C
    n = Info.TrigIdx(c,m) + Cfg.AlazarLatency + (0:R-1);
    cmd = n - Cfg.TrigDelay - Cfg.AlazarLatency;
    x = W(1,cmd+1); y = W(2,cmd+1); r = hypot(x-src(1), y-src(2));
    u = uamp*exp(-r/8).*sin(w0*n*Ts - k0*r);
    phi = 4*pi*nidx*u/lam + phStatic((c-1)*R+(1:R)) + 0.01*randn(1,R);
    fr = 2000 + 600*cos(2*pi*z0*kk/nS + phi) + 5*randn(nS,R);
    fwrite(fid, uint16(round(fr)), 'uint16');
  end
end
fclose(fid);
Res = SBM_TwinOCT_ProcessStarter(fullfile(out,'smoke.bin'), fullfile(out,'sidecar.mat'), struct('MakeFigures',true));
[XX,YY] = meshgrid(Res.x, Res.y); rr = hypot(XX-src(1), YY-src(2));
Uid = exp(-rr/8).*exp(1j*(-k0*rr - pi/2));
e = angle(Res.Phasor.*conj(Uid)); e = angle(exp(1j*(e - angle(mean(exp(1j*e(:)))))));
fprintf('surface idx median %d (expect %d); phase err RMS %.3f rad; median speed (r>2mm) %.2f m/s\n', median(Res.SurfaceIdx(:)), z0+1, rms(e(:)), median(Res.Speed_ms(rr>2)));
exportgraphics(gcf, fullfile(fileparts(here),'Validation','SBM_ProcessStarter_SyntheticBin.png'),'Resolution',130);
