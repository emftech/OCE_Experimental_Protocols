function [FramesDF, Info] = MMB_DirectionalFilterOutward(Frames, Geom, fs, f_an, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_DirectionalFilterOutward
%
% Directional (propagation-direction) filter for Multimeridian MB-mode
% data. Every meridian passes through the excitation point (needle) at its
% centre, so the wave of interest travels OUTWARD along both halves of the
% meridian: towards +s for s > 0 and towards -s for s < 0. Waves reflected
% by the phantom boundaries travel inward and corrupt the phase gradient.
%
% For each meridian the space-time map v(s,t) is transformed to the
% wavenumber-frequency domain (2D FFT, kernel exp(-j 2 pi (kappa s + f t))):
%
%   a wave cos(2 pi f t - k s) moving towards +s has its energy at
%   (kappa, f) = (-k/2pi, +f) and (+k/2pi, -f)  ->  kappa*f < 0
%
% so the +s and -s propagating fields are separated with the complementary
% quadrant masks  H+ = [kappa*f < 0],  H- = [kappa*f > 0]  (same idea as
% DirectionalFilter2.m). The outward field is then stitched:
%
%   v_out(s,t) = v+(s,t) for s > 0,   v-(s,t) for s < 0
%
% Mask smoothing. Each half meridian is only ~5 mm long (about one
% wavelength), so the spectral lobe of the wave is wide along kappa and
% crosses kappa = 0. A hard quadrant mask then distorts the phase near the
% needle gap and the ends of the meridian. The masks are therefore smoothed
% with an anisotropic Gaussian whose width along kappa is SigmaKappaRel
% times the dominant wavenumber kappa0 (measured from the data at f_an), and
% SigmaFreqBins bins along frequency. Synthetic validation (8 meridians,
% 12 mm, outward wave + 60 % inward reflection + noise) gives |bias| <= ~3 %
% at 2 and 3.5 kHz with this setting (see Processing 07-10-2026 note).
%
% Rows of reversed meridians (every 2nd meridian, Geom.s_mm decreasing) are
% sorted by s before filtering and returned in acquisition order.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   Frames : [nPos*nMer x Nt] traces in acquisition order (NaN rows allowed:
%            needle footprint / missing surface)
%   Geom   : struct from OCE_MB_MeridianGeometry (s_mm, nPosPerMer, nMer)
%   fs     : temporal sampling frequency [Hz]
%   f_an   : analysis (excitation) frequency [Hz], used to find kappa0
%
% Optional name-value parameters:
%   'SigmaKappaRel' : mask smoothing along kappa, relative to kappa0 (0.5)
%   'SigmaFreqBins' : mask smoothing along frequency in bins (2)
%   'PadFactor'     : zero padding factor in both dimensions (4)
%   'KeepSpectrumOfMeridian' : meridian whose k-f spectra are returned (1)
%   'Source_mm'     : [x0 y0] excitation point (default [0 0]). Each
%                     meridian is split at the projection of the source
%                     on it, s0 = x0 cos(theta) + y0 sin(theta), so the
%                     wave is outward from the needle even if the needle
%                     is not at the centre of the meridian pattern.
%
% Outputs:
%   FramesDF : outward-filtered traces, same size/order as Frames
%   Info     : kappa0_cyc_per_m, Sigma (bins), OutwardEnergyFraction
%              [1 x nMer] (outward energy / total, a reflection indicator),
%              and the shifted k-f magnitude spectra of one meridian
%              before/after filtering (kappa_cyc_per_mm, f_Hz, Spec_Before,
%              Spec_After)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'SigmaKappaRel', 0.5);
addParameter(p, 'SigmaFreqBins', 2);
addParameter(p, 'PadFactor', 4);
addParameter(p, 'KeepSpectrumOfMeridian', 1);
addParameter(p, 'Source_mm', [0 0]);
parse(p, varargin{:});
opts = p.Results;

nPos = Geom.nPosPerMer;
nMer = Geom.nMer;
Nt = size(Frames, 2);

Ns  = 2^nextpow2(opts.PadFactor*nPos);
NtF = 2^nextpow2(opts.PadFactor*Nt);

ds_m = Geom.ds_mm*1e-3;
kap = fftFreqAxis(Ns, 1/ds_m);         % [cycles/m]
fr  = fftFreqAxis(NtF, fs);            % [Hz]
[FR, KAP] = meshgrid(fr, kap);

%% ===================== 2D spectra of all meridians =======================

Fall = complex(zeros(Ns, NtF, nMer));
sSorted = zeros(nPos, nMer);
ordAll = zeros(nPos, nMer);
badAll = false(nPos, nMer);

for b = 1:nMer
    rows = (b-1)*nPos + (1:nPos);
    [sSorted(:,b), ordAll(:,b)] = sort(Geom.s_mm(:, b));
    X = Frames(rows(ordAll(:,b)), :);
    badAll(:,b) = any(~isfinite(X), 2);
    X(badAll(:,b), :) = 0;
    Fall(:,:,b) = fft2(X, Ns, NtF);
end

%% ===================== Dominant wavenumber at f_an =======================

[~, iF] = min(abs(fr - f_an));
specK = sum(abs(Fall(:, iF, :)), 3);      % summed over meridians
specK(abs(kap) < 1/(nPos*ds_m)) = 0;      % ignore the kappa ~ 0 (DC) bins
[~, iK] = max(specK);
kappa0 = abs(kap(iK));

dKap = abs(kap(2) - kap(1));
sigK = max(opts.SigmaKappaRel*kappa0/dKap, 1);
sigF = opts.SigmaFreqBins;

Hplus  = double(KAP.*FR < 0) + 0.5*double(KAP.*FR == 0);
Hminus = double(KAP.*FR > 0) + 0.5*double(KAP.*FR == 0);
Hplus  = ifftshift(imgaussfilt(fftshift(Hplus),  [sigK sigF], 'Padding', 'replicate'));
Hminus = ifftshift(imgaussfilt(fftshift(Hminus), [sigK sigF], 'Padding', 'replicate'));

%% ===================== Filter and stitch the outward field ===============

FramesDF = nan(size(Frames));
OutFrac = nan(1, nMer);
Info = struct();

for b = 1:nMer
    rows = (b-1)*nPos + (1:nPos);
    s = sSorted(:,b) - (opts.Source_mm(1)*Geom.ux(b) + opts.Source_mm(2)*Geom.uy(b));
    bad = badAll(:,b);
    F = Fall(:,:,b);

    Xp = real(ifft2(F .* Hplus));
    Xm = real(ifft2(F .* Hminus));
    Xp = Xp(1:nPos, 1:Nt);
    Xm = Xm(1:nPos, 1:Nt);

    Y = Xm;
    Y(s > 0, :) = Xp(s > 0, :);
    Y(s == 0, :) = 0.5*(Xp(s == 0, :) + Xm(s == 0, :));

    % Energy of the inward field (reflections) vs the outward field
    Yin = Xp;
    Yin(s > 0, :) = Xm(s > 0, :);
    eOut = sum(Y(~bad, :).^2, 'all');
    eIn  = sum(Yin(~bad, :).^2, 'all');
    OutFrac(b) = eOut/(eOut + eIn + eps);

    Y(bad, :) = NaN;
    out = nan(nPos, Nt);
    out(ordAll(:,b), :) = Y;
    FramesDF(rows, :) = out;

    if b == opts.KeepSpectrumOfMeridian
        Info.kappa_cyc_per_mm = fftshift(kap)*1e-3;
        Info.f_Hz = fftshift(fr);
        Info.Spec_Before = fftshift(abs(F));
        Yz = Y; Yz(bad, :) = 0;
        Info.Spec_After = fftshift(abs(fft2(Yz, Ns, NtF)));
        Info.Meridian = b;
        Info.Hplus = fftshift(Hplus);
    end
end

Info.kappa0_cyc_per_m = kappa0;
Info.k0_rad_per_m = 2*pi*kappa0;
Info.c_kf_mps = f_an/kappa0;             % rough speed from the k-f peak
Info.Sigma_bins = [sigK sigF];
Info.OutwardEnergyFraction = OutFrac;
Info.Options = opts;

end

%% ----------------------------------------------------------------------- %
function f = fftFreqAxis(N, fsamp)
% Signed FFT frequency axis in FFT order.
f = (0:N-1).' * fsamp/N;
f(f >= fsamp/2) = f(f >= fsamp/2) - fsamp;
end
