function Res = SBM_TwinOCT_NullCheck(binFile, sidecarFile, Opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_NullCheck
%
% Processing of an SBM acquisition WITHOUT excitation (the "null" run).
% With nothing vibrating, every frame must be identical, so this run tests
% everything except the wave itself:
%   1. File   : header vs sidecar (R, C, N_F), file size.
%   2. Raster : B-scans (depth x lateral) of 3 rows, en-face max projection,
%               surface depth map. Is the geometry right?
%   3. Galvo lag (B-frame Signal Delay): rows scanned L->R and R->L are
%               displaced in opposite directions when the trigger delay does
%               not match the real galvo lag. The relative shift d between
%               even and odd rows gives the correction
%                   Delta = -d/2  [A-lines],  new delay = old delay + Delta.
%   4. Repeatability: frame-to-frame correlation and lateral shift of the
%               en-face image (galvo repeatability / sample drift).
%   5. Phase stability: surface displacement across the N_F frames, bulk
%               (common-mode) axial drift, per-pixel displacement noise
%               [nm] vs SNR, and the noise floor at the f_exc bin (bin N_P):
%               the smallest vibration amplitude SBM can detect.
%
%   Res = SBM_TwinOCT_NullCheck('PC_SBM_..._noExc_r1.bin')            % sidecar found automatically
%   Res = SBM_TwinOCT_NullCheck(bin, 'SBM_Waveform_20261017_1530.mat')
%   Res = SBM_TwinOCT_NullCheck(bin, '', struct('Cut_Depth',[30 900]))
%
% Sidecar search when sidecarFile is empty: <bin>_Sidecar.mat, <bin>.mat,
% then SBM_Waveform_LAST.mat next to the .bin. If none exists, Info is
% rebuilt from SBM_TwinOCT_Config.m with the header's trigger delay (warns:
% valid only if the config was not changed since the acquisition).
% Figures and a summary .txt go to <bin folder>/NullCheck_<bin name>/.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2, sidecarFile = ''; end
if nargin < 3, Opt = struct(); end
Opt = setDef(Opt, 'Cut_Depth', [30 750]);          % depth bins kept (as MB pipeline)
Opt = setDef(Opt, 'lambda0_m', 1.3e-6);
Opt = setDef(Opt, 'n_index', 1.39);
Opt = setDef(Opt, 'SurfaceHalfWin', 1);            % +- depth bins around surface
Opt = setDef(Opt, 'MaxShift', 8);                  % lag search range [A-lines]
Opt = setDef(Opt, 'AxialRes_um', []);              % [] -> depth in pixels
Opt = setDef(Opt, 'OutDir', '');
Opt = setDef(Opt, 'MakeFigures', true);

% Resolve full paths first: fopen also searches the MATLAB path, dir/isfile
% and memmapfile only the current folder (cause of the first Twin error).
binFile = resolveFile(binFile, 'Raw .bin file');
[fdir, fname] = fileparts(binFile);
if isempty(Opt.OutDir), Opt.OutDir = fullfile(fdir, ['NullCheck_' fname]); end
if ~exist(Opt.OutDir, 'dir'), mkdir(Opt.OutDir); end

%% 1. File, header, sidecar --------------------------------------------------
[Hd, off] = readHeader(binFile);
[Info, scSrc] = loadInfo(binFile, sidecarFile, Hd);
nS = Hd.samples;  R = Hd.Alines;  C = Hd.Bframes;  nF = Hd.Repeats;
L = struct();
L.Header = sprintf('samples %d | A-lines/B-frame %d | B-frames %d | 3D repeats %d | delay %d | %g x %g mm', ...
    nS, R, C, nF, Hd.Delay, Hd.Hlen_mm, Hd.Vlen_mm);
L.Sidecar = scSrc;
L.HeaderMatchesSidecar = R == Info.R && C == Info.C && nF == Info.Cfg.nFrames;
if ~L.HeaderMatchesSidecar
    error('SBM:Header', ['Header (R=%d, C=%d, N_F=%d) does not match the sidecar ' ...
        '(R=%d, C=%d, N_F=%d). Wrong sidecar or GUI values?'], R, C, nF, ...
        Info.R, Info.C, Info.Cfg.nFrames);
end
if isfinite(Hd.Delay) && Hd.Delay ~= Info.Cfg.TrigDelay
    warning('SBM:Delay', 'Header delay %d ~= sidecar delay %d; using the sidecar.', ...
        Hd.Delay, Info.Cfg.TrigDelay);
end
d = dir(binFile);
L.ExpectedBytes = off + 2*nS*R*C*nF;  L.ActualBytes = d.bytes;
L.FileComplete = d.bytes >= L.ExpectedBytes;
if ~L.FileComplete
    error('SBM:Truncated', 'File has %d bytes, %d expected (acquisition interrupted?).', ...
        d.bytes, L.ExpectedBytes);
end
L.ExcModeInSidecar = string(Info.Cfg.ExcMode);
fprintf('\n=== SBM null check: %s ===\n%s\nSidecar: %s\n', fname, L.Header, scSrc);

m = memmapfile(binFile, 'Offset', off, 'Format', {'uint16', [nS R*C nF], 'raw'});
Nper = R*C;  nA = Info.Cfg.nAscans;  nB = Info.Cfg.nBscans;
PI = Info.PixelIdx;                                     % (b,col) -> A-line in frame record
img = @(v) reshape(v(PI), nB, nA);

%% 2. Complex signal, structure, surface ------------------------------------
probe = unique(round(linspace(1, nF, min(nF, 8))));
DC = zeros(nS, numel(probe));
for i = 1:numel(probe)
    DC(:,i) = median(double(m.Data.raw(:,:,probe(i))), 2);
end
DC = median(DC, 2);
hw = hann(nS);
zIdx = Opt.Cut_Depth(1):min(Opt.Cut_Depth(2), nS/2);
nZ = numel(zIdx);

Amp = zeros(nZ, Nper);
for f = 1:nF
    Cx = fft(hw .* (double(m.Data.raw(:,:,f)) - DC));
    Amp = Amp + abs(Cx(zIdx,:));
end
Amp = Amp / nF;
noiseAmp = median(Amp(max(1,end-49):end, :), 'all');    % deepest 50 bins
[~, zs] = max(movmean(Amp, 5, 1), [], 1);
zs = round(medfilt1(zs, 5));

w = -Opt.SurfaceHalfWin:Opt.SurfaceHalfWin;
Ssurf = complex(zeros(Nper, nF, numel(w)));
for f = 1:nF
    Cx = fft(hw .* (double(m.Data.raw(:,:,f)) - DC));
    Cx = Cx(zIdx,:);
    for iw = 1:numel(w)
        zz = min(max(zs + w(iw), 1), nZ);
        Ssurf(:,f,iw) = Cx(sub2ind(size(Cx), zz, 1:Nper)).';
    end
end
SurfAmpF = max(abs(Ssurf), [], 3);                      % [Nper x nF]
SurfAmp  = mean(SurfAmpF, 2);
SNRdB    = 20*log10(SurfAmp / noiseAmp);

% axes
x_mm = linspace(-Hd.Hlen_mm/2, Hd.Hlen_mm/2, nA);
y_mm = linspace( Hd.Vlen_mm/2,-Hd.Vlen_mm/2, nB);
if isempty(Opt.AxialRes_um)
    zAx = zIdx;  zLab = 'Depth [px]';  zsc = 1;
else
    zsc = Opt.AxialRes_um*1e-3/Opt.n_index;  zAx = zIdx*zsc;  zLab = 'Depth [mm]';
end
SurfDepth = img(zs + zIdx(1) - 1) * zsc;
MIP = img(20*log10(max(Amp, [], 1)));
SurfImg = img(20*log10(SurfAmp));

%% 3. Galvo lag from the serpentine (odd vs even rows) -----------------------
Dg = Info.Cfg.TrigDelay;
cand = {'surface amplitude', SurfImg; 'max projection', MIP; 'surface depth', SurfDepth};
Lag = struct('Image', {}, 'd', {}, 'Peak', {}, 'Curve', {});
for c = 1:size(cand,1)
    [dd, pk, curve, sh] = rowShift(cand{c,2}, Opt.MaxShift);
    Lag(c) = struct('Image', cand{c,1}, 'd', dd, 'Peak', pk, 'Curve', curve);
end
% combined score (mean over images): depth disambiguates periodic reflectivity
combo = mean(vertcat(Lag.Curve), 1);
[pkC, iC] = max(combo);
dBest = sh(iC);
[~, best] = max([Lag.Peak]);                            % image shown in the figure
Delta = -dBest/2;
L.RowShift_Alines = dBest;
L.LagCorrection_Alines = Delta;
L.RecommendedDelay = Dg + round(Delta);
L.LagEstimateFrom = 'combined (surface amplitude + max projection + surface depth)';
L.LagScorePeak = pkC;

%% 4. Repeatability across frames --------------------------------------------
ref = SurfImg;
rho = zeros(nF,1);  shiftXY = zeros(nF,2);
for f = 1:nF
    A = img(20*log10(SurfAmpF(:,f)));
    rho(f) = corr2d(A, ref);
    shiftXY(f,:) = xcorrShift(A, ref);
end
L.MinFrameCorr = min(rho);
L.MaxFrameShift_px = max(abs(shiftXY), [], 'all');

%% 5. Phase stability ----------------------------------------------------------
ph = angle(sum(Ssurf .* conj(mean(Ssurf, 2)), 3));
ph = unwrap(ph, [], 2);
k2nm = Opt.lambda0_m/(4*pi*Opt.n_index)*1e9;            % nm per rad
valid = PI(:);
tTrue = (Info.FrameStart - Info.FrameStart(1)) * Info.Cfg.Ts_OCT;     % s
bulkRad = median(ph(valid,:), 1);                       % common mode per frame (rad)
bulk_nm = (bulkRad - bulkRad(1))' * k2nm;               % axial drift / table vibration
tq = tTrue / max(tTrue(end), eps);
Bm = [ones(nF,1), tq];
fitRes = @(P) (P - (Bm*(Bm \ P.')).') * k2nm;          % per-pixel offset + linear removed
rawRes_nm = fitRes(ph);                                 % without bulk correction
res_nm = fitRes(ph - bulkRad);                          % after bulk-motion correction
sig_nm = std(res_nm, 0, 2);
sigRaw_nm = std(rawRes_nm, 0, 2);
NP = Info.Cfg.nPeriods;
Uk = (2/nF) * abs(fft(res_nm(valid,:), [], 2));          % [Npix x nF] amplitude spectrum
specMean = mean(Uk(:, 1:floor(nF/2)+1), 1);
L.DisplNoise_nm_median = median(sig_nm(valid));           % after bulk correction
L.DisplNoiseNoBulkCorr_nm = median(sigRaw_nm(valid));
L.FloorAt_fexc_nm = mean(Uk(:, NP+1));
L.BulkDrift_nm_pp = max(bulk_nm) - min(bulk_nm);
L.SurfaceSNR_dB_median = median(SNRdB(valid));
L.SurfaceDepth_median = median(SurfDepth(:));

fprintf(['File OK (%.1f MB) | surface SNR %.1f dB | displacement noise %.2f nm (median, ' ...
    'bulk-corrected; %.2f nm without) | floor at f_exc bin %.3f nm | bulk drift %.1f nm p-p\n'], ...
    L.ActualBytes/2^20, L.SurfaceSNR_dB_median, L.DisplNoise_nm_median, ...
    L.DisplNoiseNoBulkCorr_nm, L.FloorAt_fexc_nm, L.BulkDrift_nm_pp);
fprintf('Frame repeatability: min corr %.3f, max lateral shift %.2f px\n', ...
    L.MinFrameCorr, L.MaxFrameShift_px);
fprintf(['Serpentine shift (even - odd rows) = %+.2f A-lines (combined score %.2f; per image: %s)\n' ...
    '  -> galvo lag = delay %+.2f A-lines -> set B-frame Signal Delay = %d (now %d)\n'], ...
    dBest, pkC, mat2str([Lag.d], 3), Delta, L.RecommendedDelay, Dg);

Res = struct('Summary', L, 'Header', Hd, 'Info', Info, 'x_mm', x_mm, 'y_mm', y_mm, ...
    'SurfDepth', SurfDepth, 'SurfImg_dB', SurfImg, 'MIP_dB', MIP, 'SNR_dB', img(SNRdB), ...
    'DisplNoise_nm', img(sig_nm), 'Bulk_nm', bulk_nm, 't_s', tTrue, 'FrameCorr', rho, ...
    'FrameShift_px', shiftXY, 'Lag', Lag, 'Spectrum_nm', specMean, ...
    'Residual_nm', res_nm, 'AmpMean', Amp, 'zIdx', zIdx);

%% Figures ---------------------------------------------------------------------
if Opt.MakeFigures
    Ad = 20*log10(Amp);  cl = prctile(Ad(:), [5 99.7]);
    f1 = figure('Name','Null check: structure','Color','w','Position',[40 60 1350 720]);
    t1 = tiledlayout(f1, 2, 3, 'TileSpacing','compact');
    rows = unique([1, ceil(nB/2), nB]);
    for r = rows
        nexttile(t1);
        Bsc = Ad(:, PI(r,:));                           % depth x columns (spatial order)
        imagesc(x_mm, zAx, Bsc); clim(cl); colormap(gca, gray);
        hold on; plot(x_mm, SurfDepth(r,:), 'r.');
        xlabel('x [mm]'); ylabel(zLab); title(sprintf('B-scan row %d (y = %.2f mm)', r, y_mm(r)));
    end
    nexttile(t1); imagesc(x_mm, y_mm, MIP); axis image; set(gca,'YDir','normal');
    colormap(gca, gray); colorbar; title('En-face max projection [dB]'); xlabel('x [mm]'); ylabel('y [mm]');
    nexttile(t1); imagesc(x_mm, y_mm, SurfDepth); axis image; set(gca,'YDir','normal');
    colormap(gca, parula); colorbar; title(['Surface ' lower(zLab)]);
    xlabel('x [mm]'); ylabel('y [mm]');
    nexttile(t1); imagesc(x_mm, y_mm, img(SNRdB)); axis image; set(gca,'YDir','normal');
    colormap(gca, turbo); colorbar; title('Surface SNR [dB]'); xlabel('x [mm]'); ylabel('y [mm]');
    title(t1, sprintf('%s | %d x %d A-lines, %d frames', fname, nA, nB, nF), 'Interpreter','none');
    saveFig(f1, Opt.OutDir, 'NullCheck_1_Structure');

    f2 = figure('Name','Null check: galvo lag','Color','w','Position',[60 80 1350 650]);
    t2 = tiledlayout(f2, 2, 3, 'TileSpacing','compact');
    I0 = cand{best,2};
    nexttile(t2); imagesc(x_mm, y_mm, I0); axis image; set(gca,'YDir','normal'); colorbar;
    title(sprintf('As acquired (%s)', cand{best,1})); xlabel('x [mm]'); ylabel('y [mm]');
    nexttile(t2); imagesc(x_mm, y_mm, shiftEven(I0, -dBest)); axis image; set(gca,'YDir','normal'); colorbar;
    title(sprintf('Even rows shifted by %+.2f A-lines (preview of the fix)', -dBest));
    xlabel('x [mm]'); ylabel('y [mm]');
    nexttile(t2); hold on;
    for c = 1:numel(Lag), plot(sh, Lag(c).Curve, 'LineWidth', 1.1); end
    plot(sh, combo, 'k-', 'LineWidth', 2.2);
    xline(dBest, 'k--'); grid on; xlabel('Shift of even rows vs odd rows [A-lines]');
    ylabel('Match score  1 - MSD/var'); legend([{Lag.Image}, {'combined'}], 'Location','best');
    title(sprintf('d = %+.2f -> lag correction %+.2f -> delay %d', dBest, Delta, L.RecommendedDelay));
    nexttile(t2, [1 3]);
    rec = 1:min(Nper, 3*Info.S + Info.Cfg.nAscans);       % first 3 B-scan steps of the record
    isA = false(Nper,1); isA(PI(:)) = true;
    aR = 20*log10(SurfAmp(rec));  vR = rec(isA(rec));
    yyaxis left; plot(rec, aR, 'k-'); hold on; plot(vR, aR(isA(rec)), 'g.', 'MarkerSize', 9);
    ylabel('Surface amplitude [dB]');
    yyaxis right; plot(rec, zs(rec) + zIdx(1) - 1, 'b-'); ylabel('Surface depth [px]');
    xlabel('A-line index in the record'); grid on;
    legend('amplitude', 'valid A-lines', 'surface depth', 'Location', 'best');
    title('Start of the record (green = valid A-lines, the rest are turn samples): with the right delay, features are mirror-symmetric about each turn');
    saveFig(f2, Opt.OutDir, 'NullCheck_2_GalvoLag');

    f3 = figure('Name','Null check: repeatability and phase stability','Color','w','Position',[80 100 1350 720]);
    t3 = tiledlayout(f3, 2, 3, 'TileSpacing','compact');
    nexttile(t3); plot(1:nF, rho, 'k.-'); grid on; xlabel('Frame'); ylabel('Correlation with mean');
    title('En-face repeatability'); ylim([min(0.9, min(rho)) 1]);
    nexttile(t3); plot(1:nF, shiftXY(:,1), 'b.-', 1:nF, shiftXY(:,2), 'r.-'); grid on;
    xlabel('Frame'); ylabel('Shift [px]'); legend('x','y'); title('Lateral frame shift (galvo / drift)');
    nexttile(t3); plot(tTrue*1e3, bulk_nm, 'k.-'); grid on; xlabel('Time [ms]'); ylabel('nm');
    title(sprintf('Bulk axial displacement (%.1f nm p-p)', L.BulkDrift_nm_pp));
    nexttile(t3); imagesc(x_mm, y_mm, img(sig_nm)); axis image; set(gca,'YDir','normal');
    colormap(gca, turbo); colorbar; clim([0 prctile(sig_nm(valid), 98)]);
    title(sprintf('Displacement noise [nm] (bulk-corrected), median %.2f', L.DisplNoise_nm_median));
    xlabel('x [mm]'); ylabel('y [mm]');
    nexttile(t3);
    snr = SNRdB(valid);  sg = sig_nm(valid);
    semilogy(snr, sg, '.', 'Color', [.5 .5 .5]); hold on;
    sv = linspace(min(snr), max(snr), 100);
    sigma_n = noiseAmp/sqrt(2*log(2));                  % Rayleigh median -> per-quadrature sigma
    semilogy(sv, k2nm*sigma_n./(noiseAmp*10.^(sv/20)), 'r-', 'LineWidth', 1.5);
    grid on; xlabel('Surface SNR [dB]'); ylabel('\sigma_u [nm]');
    legend('pixels', 'single-bin SNR limit', 'Location','best'); title('Noise vs SNR');
    nexttile(t3); stem(0:numel(specMean)-1, specMean, 'k', 'filled', 'MarkerSize', 3); hold on;
    stem(NP, specMean(NP+1), 'r', 'filled'); grid on;
    xlabel('Frame-DFT bin (bin N_P = f_{exc})'); ylabel('Mean amplitude [nm]');
    title(sprintf('Noise floor at f_{exc}: %.3f nm (no excitation)', L.FloorAt_fexc_nm));
    saveFig(f3, Opt.OutDir, 'NullCheck_3_Stability');
end

save(fullfile(Opt.OutDir, 'NullCheck_Results.mat'), 'Res', '-v7.3');
writeSummary(fullfile(Opt.OutDir, 'NullCheck_Summary.txt'), fname, L);
fprintf('Results in %s\n', Opt.OutDir);
end

%% ==========================================================================
function [d, pk, curve, sh] = rowShift(I, maxS)
% Relative x-shift of even rows vs the mean of their odd neighbours.
% Score = 1 - MSD/var(ref): unlike a correlation coefficient it is sensitive
% to offsets, so a smoothly tilted surface still has a sharp optimum.
[nB, nA] = size(I);
sh = -maxS:0.1:maxS;  curve = zeros(size(sh));  nUse = 0;
v0 = var(I(isfinite(I)));
for b = 2:2:nB-1
    ref = mean(I([b-1 b+1], :), 1);
    pr  = I(b, :);
    if any(~isfinite([ref pr])), continue; end
    nUse = nUse + 1;
    for k = 1:numel(sh)
        q = interp1(1:nA, pr, (1:nA) + sh(k), 'linear', NaN);
        ok = isfinite(q);
        if nnz(ok) < 5, curve(k) = curve(k) - 1; continue; end
        curve(k) = curve(k) + 1 - mean((q(ok) - ref(ok)).^2)/(v0 + eps);
    end
end
curve = curve / max(nUse, 1);
[pk, i] = max(curve);
d = sh(i);          % even row content sits at column (ref col + d)
end

function J = shiftEven(I, s)
J = I;  nA = size(I,2);
for b = 2:2:size(I,1)
    J(b,:) = interp1(1:nA, I(b,:), (1:nA) - s, 'linear', NaN);
end
end

function r = corr2d(A, B)
a = A(:) - mean(A(:));  b = B(:) - mean(B(:));
r = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2) + eps);
end

function s = xcorrShift(A, B)
% sub-pixel [dx dy] of A relative to B (FFT cross-correlation + parabola)
A = A - mean(A(:));  B = B - mean(B(:));
sz = 2*size(A);
Xc = real(ifft2(fft2(A, sz(1), sz(2)) .* conj(fft2(B, sz(1), sz(2)))));
Xc = fftshift(Xc);
[~, i] = max(Xc(:));  [iy, ix] = ind2sub(size(Xc), i);
iy = min(max(iy, 2), size(Xc,1)-1);  ix = min(max(ix, 2), size(Xc,2)-1);
c0 = floor(sz/2) + 1;
s = [ix - c0(2) + par(Xc(iy, ix-1:ix+1)), iy - c0(1) + par(Xc(iy-1:iy+1, ix))];
    function o = par(v)
        if numel(v) < 3 || any(~isfinite(v)), o = 0; return; end
        den = v(1) - 2*v(2) + v(3);
        o = 0; if abs(den) > eps, o = 0.5*(v(1) - v(3))/den; end
    end
end

function [Info, src] = loadInfo(binFile, sidecarFile, Hd)
[fdir, fname] = fileparts(binFile);
if ~isempty(sidecarFile)
    % explicit sidecar: must exist (as given, on the path, or next to the .bin)
    src = resolveFile(sidecarFile, 'Sidecar .mat', fdir);
    Info = readInfo(src);
    return
end
cands = {fullfile(fdir, [fname '_Sidecar.mat']), fullfile(fdir, [fname '.mat']), ...
    fullfile(fdir, 'SBM_Waveform_LAST.mat')};
for i = 1:numel(cands)
    if isfile(cands{i})
        Info = readInfo(cands{i});  src = cands{i};  return
    end
end
warning('SBM:NoSidecar', ['No sidecar found next to the .bin (looked for %s_Sidecar.mat, ' ...
    '%s.mat, SBM_Waveform_LAST.mat in %s): rebuilding Info from SBM_TwinOCT_Config.m ' ...
    '(valid only if the config is the one used for this acquisition).'], fname, fname, fdir);
Cfg = SBM_TwinOCT_Config();
if isfinite(Hd.Delay), Cfg.TrigDelay = Hd.Delay; end
[~, Info] = SBM_TwinOCT_Waveform(Cfg);
src = 'rebuilt from SBM_TwinOCT_Config.m';
end

function Info = readInfo(file)
vars = who('-file', file);
if ~ismember('Info', vars)
    error('SBM:Sidecar', ['%s has no variable "Info" (it has: %s). Use the ' ...
        'SBM_Waveform_*.mat written by the LabVIEW script node.'], file, strjoin(vars, ', '));
end
S = load(file, 'Info');  Info = S.Info;
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

function [Hd, off] = readHeader(binFile)
fid = fopen(binFile, 'r');
if fid < 0, error('SBM:File', 'Cannot open %s', binFile); end
Lh = fread(fid, 1, 'uint16');
fseek(fid, 4, 'bof');
H = fread(fid, Lh, '*char')';
fclose(fid);
off = 4 + Lh;
tok = regexp(H, '([^:;\r\n]+):\s*([^;\r\n]*)', 'tokens');
key = strtrim(cellfun(@(t) t{1}, tok, 'UniformOutput', false));
val = cellfun(@(t) strtrim(t{2}), tok, 'UniformOutput', false);
    function v = get(k)
        i = find(contains(key, k), 1);
        if isempty(i), v = NaN; else, v = str2double(val{i}); end
    end
Hd = struct('samples', get('spectral samples per A-line'), ...
    'Alines', get('A-lines per B-frame'), 'Bframes', get('B-frames in 3D scan'), ...
    'Repeats', get('3D repeated scans'), 'Delay', get('B-frame trigger delay'), ...
    'Hlen_mm', get('Horizontal scan length'), 'Vlen_mm', get('Vertical scan length'), ...
    'Raw', H);
end

function saveFig(f, d, name)
exportgraphics(f, fullfile(d, [name '.png']), 'Resolution', 130);
savefig(f, fullfile(d, [name '.fig']));
end

function writeSummary(file, fname, L)
fid = fopen(file, 'w');
fprintf(fid, 'SBM null check: %s\n%s\n\n', fname, char(datetime('now')));
fn = fieldnames(L);
for i = 1:numel(fn)
    v = L.(fn{i});
    if isnumeric(v) || islogical(v), v = mat2str(v, 5); end
    fprintf(fid, '%-24s : %s\n', fn{i}, string(v));
end
fclose(fid);
end

function S = setDef(S, f, v)
if ~isfield(S, f) || isempty(S.(f)), S.(f) = v; end
end
