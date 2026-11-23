function Res = SBM_TwinOCT_ProcessStarter(binFile, sidecarFile, Opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_ProcessStarter
%
% First-pass processing of one SBM acquisition (raw .bin + sidecar .mat
% saved by the LabVIEW node), meant to check the protocol on the bench:
%   1. header + size check (layout R x C x N_F A-lines, MB-style),
%   2. complex OCT signal (DC-spectrum subtraction, Hann, FFT; same math as
%      ReadCplx_Volume_BS / MMB_ComputeCplx), one frame at a time,
%   3. surface depth per raster position (max of the frame-averaged
%      amplitude), complex surface signal per A-line,
%   4. surface displacement per frame  u_m = lambda0/(4 pi n) *
%      angle(s_m conj(<s>_m)), linear drift removed,
%   5. SBM_TwinOCT_Reorganize -> SBM_TwinOCT_Recover (Eq. 59),
%   6. quick phase-gradient speed map (central differences of the
%      recovered phasor). Use the thesis 2D PG / PG-TV estimators
%      (Prepare_SBM_2DPhaseGradient + PhaseGradient2D_TV_FromSims) or
%      OCE_PG_* for the real elastograms.
%
%   Res = SBM_TwinOCT_ProcessStarter('PC_SBM_....bin', 'SBM_Waveform_....mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3, Opt = struct(); end
Opt = setDef(Opt, 'Cut_Depth', [30 750]);
Opt = setDef(Opt, 'lambda0_m', 1.3e-6);
Opt = setDef(Opt, 'n_index', 1.39);
Opt = setDef(Opt, 'SurfaceHalfWin', 1);
Opt = setDef(Opt, 'FOV_mm', []);                 % [Lx Ly] for mm axes
Opt = setDef(Opt, 'MakeFigures', true);

% full paths (fopen searches the MATLAB path; dir and memmapfile do not)
binFile = resolveFile(binFile, 'Raw .bin file');
sidecarFile = resolveFile(sidecarFile, 'Sidecar .mat', fileparts(binFile));
vars = who('-file', sidecarFile);
if ~ismember('Info', vars)
    error('SBM:Sidecar', '%s has no variable "Info" (it has: %s).', sidecarFile, strjoin(vars, ', '));
end
S = load(sidecarFile, 'Info');  Info = S.Info;
[Hd, off] = readHeader(binFile);
nS = Hd.samples;  R = Hd.Alines;  C = Hd.Bframes;  nF = Hd.Repeats;
if R ~= Info.R || C ~= Info.C || nF ~= Info.Cfg.nFrames
    error('SBM:Header', ['Header (R=%d, C=%d, N_F=%d) does not match the sidecar ' ...
        '(R=%d, C=%d, N_F=%d).'], R, C, nF, Info.R, Info.C, Info.Cfg.nFrames);
end
d = dir(binFile);
need = off + 2*nS*R*C*nF;
if d.bytes < need
    error('SBM:Truncated', 'File has %d bytes, %d expected.', d.bytes, need);
end
m = memmapfile(binFile, 'Offset', off, 'Format', {'uint16', [nS R*C nF], 'raw'});

%% DC spectrum and complex signal ------------------------------------------
probe = unique(round(linspace(1, nF, min(nF, 8))));
DC = zeros(nS, numel(probe));
for i = 1:numel(probe)
    DC(:,i) = median(double(m.Data.raw(:,:,probe(i))), 2);
end
DC = median(DC, 2);
hw = hann(nS);
zIdx = Opt.Cut_Depth(1):min(Opt.Cut_Depth(2), nS/2);
Nper = R*C;

Amp = zeros(numel(zIdx), Nper);
for f = 1:nF
    Cx = fft(hw .* (double(m.Data.raw(:,:,f)) - DC));
    Amp = Amp + abs(Cx(zIdx,:));
end
Amp = Amp / nF;
[~, zs] = max(movmean(Amp, 5, 1), [], 1);                % surface per A-line position
zs = round(medfilt1(zs, 5));

% Depth bins around the surface are kept separately: adjacent bins of a
% Hann-windowed sharp reflector are anti-phase, so summing them cancels.
w = -Opt.SurfaceHalfWin:Opt.SurfaceHalfWin;
Ssurf = complex(zeros(Nper, nF, numel(w)));
for f = 1:nF
    Cx = fft(hw .* (double(m.Data.raw(:,:,f)) - DC));
    Cx = Cx(zIdx,:);
    for iw = 1:numel(w)
        zz = min(max(zs + w(iw), 1), numel(zIdx));
        Ssurf(:,f,iw) = Cx(sub2ind(size(Cx), zz, 1:Nper)).';
    end
end

%% Displacement per frame ---------------------------------------------------
% amplitude-weighted mean of the per-depth phase change w.r.t. the
% frame-averaged (static) signal
ph = angle(sum(Ssurf .* conj(mean(Ssurf, 2)), 3));
ph = unwrap(ph, [], 2);
% Slow drift between frames: fit [1, t, cos, sin] jointly and remove only
% the offset + ramp (a plain detrend would also eat part of the harmonic,
% because a sine over N_P periods is not orthogonal to a ramp).
tq = (0:nF-1)'/nF;
Bm = [ones(nF,1), tq, cos(2*pi*Info.Cfg.nPeriods*tq), sin(2*pi*Info.Cfg.nPeriods*tq)];
coef = Bm \ ph.';
ph = ph - (Bm(:,1:2)*coef(1:2,:)).';
u_nm = ph * Opt.lambda0_m/(4*pi*Opt.n_index) * 1e9;

[Frames, Ax] = SBM_TwinOCT_Reorganize(reshape(u_nm, 1, []), Info);
[Rec, U, PhiAcq, Uraw] = SBM_TwinOCT_Recover(Frames, Info);

%% Quick speed map -----------------------------------------------------------
if isempty(Opt.FOV_mm)
    if ~isempty(Info.Cfg.HVC)
        x = Ax.x_mm;  y = Ax.y_mm;
    else
        x = Ax.x_V;   y = Ax.y_V;   warning('SBM:Units', 'No HVC/VVC or FOV_mm: speed in V-units.');
    end
else
    x = linspace(-Opt.FOV_mm(1)/2, Opt.FOV_mm(1)/2, numel(Ax.x_V));
    y = linspace( Opt.FOV_mm(2)/2,-Opt.FOV_mm(2)/2, numel(Ax.y_V));
end
dx = mean(diff(x));  dy = mean(abs(diff(y)));
kx = nan(size(U));  ky = kx;                            % borders undefined
kx(:,2:end-1) = angle(U(:,3:end).*conj(U(:,1:end-2)))/(2*dx);
ky(2:end-1,:) = angle(U(1:end-2,:).*conj(U(3:end,:)))/(2*dy);
kmag = movmedian(movmedian(hypot(kx, ky), 3, 1, 'omitnan'), 3, 2, 'omitnan');
c_ms = 2*pi*Info.f_exc ./ (kmag*1e3);

Res = struct('Header', Hd, 'Info', Info, 'SurfaceIdx', zs + zIdx(1) - 1, ...
    'Frames_nm', Frames, 'Recovered_nm', Rec, 'Phasor', U, 'PhasorRaw', Uraw, ...
    'PhiAcq', PhiAcq, 'Speed_ms', c_ms, 'x', x, 'y', y, 't_equiv_s', Ax.t_equiv_s);

if Opt.MakeFigures
    figure('Color','w','Position',[80 80 1300 380]);
    tl = tiledlayout(1,4,'TileSpacing','compact');
    sA = Amp(sub2ind(size(Amp), zs, 1:Nper));
    nexttile(tl); imagesc(x, y, sA(Info.PixelIdx));
    axis image; set(gca,'YDir','normal'); title('Surface amplitude'); colorbar; colormap(gca, gray);
    nexttile(tl); imagesc(x, y, angle(Uraw)); axis image; set(gca,'YDir','normal');
    clim([-pi pi]); title('Raw phase at f_{exc}'); colormap(gca, SBM_PhaseCmap);
    nexttile(tl); imagesc(x, y, angle(U)); axis image; set(gca,'YDir','normal');
    clim([-pi pi]); title('Recovered phase (Eq. 59)'); colormap(gca, SBM_PhaseCmap);
    nexttile(tl); imagesc(x, y, c_ms); axis image; set(gca,'YDir','normal');
    clim(prctile(c_ms(:), [5 95])); title('Quick PG speed [m/s]'); colorbar; colormap(gca, turbo);
end
end

%% ==========================================================================
function [Hd, off] = readHeader(binFile)
fid = fopen(binFile, 'r');
L = fread(fid, 1, 'uint16');
fseek(fid, 4, 'bof');
H = fread(fid, L, '*char')';
fclose(fid);
off = 4 + L;
tok = regexp(H, '([^:;\r\n]+):\s*([^;\r\n]*)', 'tokens');
key = strtrim(cellfun(@(t) t{1}, tok, 'UniformOutput', false));
val = cellfun(@(t) strtrim(t{2}), tok, 'UniformOutput', false);
get = @(k) str2double(val{find(contains(key, k), 1)});
Hd = struct('samples', get('spectral samples per A-line'), ...
    'Alines', get('A-lines per B-frame'), 'Bframes', get('B-frames in 3D scan'), ...
    'Repeats', get('3D repeated scans'), 'Raw', H);
end

function full = resolveFile(f, what, altDir)
% Full path of f: as given (relative to the current folder), on the MATLAB
% path, or (name only) inside altDir.
f = char(f);  cand = {f};
if nargin > 2 && ~isempty(altDir)
    [~, n, e] = fileparts(f);  cand{end+1} = fullfile(altDir, [n e]);
end
for i = 1:numel(cand)
    if isfile(cand{i})
        d = dir(cand{i});  full = fullfile(d(1).folder, d(1).name);  return
    end
    w = which(cand{i});
    if ~isempty(w) && isfile(w), full = w;  return; end
end
error('SBM:NotFound', ['%s not found: %s\nCurrent folder: %s\n' ...
    'Pass the full path, or cd to the data folder.'], what, strjoin(cand, ' | '), pwd);
end

function S = setDef(S, f, v)
if ~isfield(S, f) || isempty(S.(f)), S.(f) = v; end
end
