function H = SBM_HarmonicRecovery(Frames, Info, kind, BandOpts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_HarmonicRecovery
%
% Harmonic isolation and acquisition-phase compensation of SBM frames
% (thesis Sec. 5.2, Eqs. 56-59), with the EXACT acquisition instant of every
% pixel taken from the hardware timeline (sidecar Info), not from a nominal
% raster model.
%
% For pixel (b,a) of frame m the A-line is recorded at the global sample
%   n = s_1 + (m-1) F + tau(b,a),   F = K P + N_phi,
% so the frame series is the harmonic sampled at the equivalent times
% (m-1) T_SBM, delayed by phi_acq(b,a) = omega (s_1 + tau(b,a)) Ts (Eq. 57).
% N_F frames span exactly N_P periods, so f_exc is EXACTLY DFT bin N_P
% (the narrowest possible band-pass, no window, no leakage):
%   U_raw = (2/N_F) sum_m v_m exp(-j 2 pi N_P (m-1)/N_F)          (Eq. 58)
%   U     = U_raw exp(-j phi_acq)                                  (Eq. 59)
% For particle velocity from a lag-1 Loupas estimate (MotionType 2) the
% sample m holds the phase change between frames m and m+1, i.e. it is
% centred half a step later: that constant offset is also compensated.
%
% Inputs:
%   Frames : [N_B x N_A x N_F] (SBM_ToFrames, y ascending)
%   Info   : sidecar Info
%   kind   : 'velocity' (lag-1, default) or 'displacement'
%   BandOpts : temporal band selection (SBM_SelectBand), as the MB
%              narrowband step. Mode 'exact' (default): bin N_P only.
%              Mode 'peak': highest spectral peak near f_exc (search window
%              f_exc*(1 -/+ SearchRel)), zero-phase Gaussian band-pass of
%              its -3 dB width (capped at MaxBW_Hz) applied circularly along
%              the frames, phasor and phi_acq evaluated at the analysis
%              frequency f_an (peak or f_exc, UsePeakFrequency).
% Output struct H:
%   PhasorRaw, Phasor : [N_B x N_A] complex (before / after compensation)
%   PhiAcq            : [N_B x N_A] accumulated acquisition phase [rad] at f_an
%   Band, f_an, BinAn : selected band (SBM_SelectBand), analysis frequency
%                       and its (possibly fractional) bin
%   Filtered          : band-passed frames [N_B x N_A x N_F]
%   Recovered         : [N_B x N_A x N_F] harmonic movie after compensation
%   Spectrum          : [1 x N_F/2+1] mean |DFT| over pixels (bins 0..N_F/2)
%   f_equiv_Hz        : equivalent frequency of every bin (bin N_P = f_exc)
%   Detect            : spatial-coherence detection test (CohExc, CohNull,
%                       Threshold = P95 of the null bins, Detected)
%   BinExc, SNR_dB    : excitation bin and its level over the median of the
%                       other non-DC bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3 || isempty(kind), kind = 'velocity'; end
if nargin < 4 || isempty(BandOpts), BandOpts = struct('Mode', 'exact'); end
nF = size(Frames, 3);
NP = Info.Cfg.nPeriods;
df = Info.f_exc / NP;                           % bin spacing = equivalent-time resolution [Hz]
m  = reshape(0:nF-1, 1, 1, []);
X  = double(Frames);
X  = X - mean(X, 3, 'omitnan');
bad = any(~isfinite(X), 3);                     % invalid pixels stay NaN in the phasor
X(~isfinite(X)) = 0;
H = struct();

% ---- temporal band (MMB_SelectNarrowBand rule, in equivalent time)
H.Band = SBM_SelectBand(X(repmat(~bad, 1, 1, nF)), nnz(~bad), Info.f_exc, NP, BandOpts);
beta = H.Band.f_an / df;                        % analysis bin (fractional if the peak is)
H.f_an = H.Band.f_an;
H.BinAn = beta;
if H.Band.Filtered
    % zero-phase Gaussian band-pass (MMB_GaussianBandpass), CIRCULAR along
    % the frames: the N_F frames span exactly N_P excitation periods
    k = (0:nF-1)';  k(k > nF/2) = k(k > nF/2) - nF;
    f = reshape(abs(k) * df, 1, 1, []);
    sigma = (H.Band.BW_Hz/2) / sqrt(log(2));
    X = real(ifft(fft(X, [], 3) .* exp(-(f - H.f_an).^2 / (2*sigma^2)), [], 3));
end

% ---- phasor at the analysis frequency and Eq. 59
PhiExc = flipud(Info.PhiAcq);                   % y ascending, as SBM_ToFrames; omega_exc (s1+tau) Ts
H.PhasorRaw = (2/nF) * sum(X .* exp(-1j*2*pi*beta*m/nF), 3);
H.PhiAcq = PhiExc * beta/NP;                    % omega_an (s1+tau) Ts
if strcmpi(kind, 'velocity')
    H.PhiAcq = H.PhiAcq + pi*beta/nF;           % lag-1 sample centred half a step later
end
H.PhasorRaw(bad) = NaN;
H.Phasor = H.PhasorRaw .* exp(-1j*H.PhiAcq);
H.Recovered = real(H.Phasor .* exp(1j*2*pi*beta*m/nF));
H.Filtered = X;

Fall = abs(fft(double(Frames) - mean(double(Frames), 3, 'omitnan'), [], 3)) * 2/nF;
nb = floor(nF/2) + 1;
H.Spectrum = squeeze(mean(reshape(Fall(:, :, 1:nb), [], nb), 1, 'omitnan'))';
H.f_equiv_Hz = (0:nb-1) * df;
H.BinExc = NP;
other = setdiff(2:nb, NP+1);
H.SNR_dB = 20*log10(H.Spectrum(NP+1) / median(H.Spectrum(other)));
[~, ipk] = max(H.Spectrum(2:nb));
H.PeakBin = ipk;                                % 0-based bin of the largest non-DC peak
H.PeakIsExcitation = (ipk == NP);

% Detection test: a real excitation-locked wave is spatially smooth after
% the phi_acq compensation (lambda ~ 3 mm >> pixel pitch), whereas noise or
% non-locked vibrations are not. Neighbour coherence of the compensated
% analysis phasor vs the same statistic on the other bins (null, unfiltered
% frames), each compensated with phi_acq*b/N_P (both signs, the larger kept).
cohN = @(U) mean([nbCoh(U, 1), nbCoh(U, 2)]);
X0 = double(Frames) - mean(double(Frames), 3, 'omitnan');  X0(~isfinite(X0)) = 0;
H.Detect.CohExc = cohN(H.Phasor);
nullBins = setdiff(1:nb-1, [floor(beta) ceil(beta) round(2*beta)]);
H.Detect.CohNull = zeros(size(nullBins));
for i = 1:numel(nullBins)
    b = nullBins(i);
    Ub = (2/nF) * sum(X0 .* exp(-1j*2*pi*b*m/nF), 3);
    H.Detect.CohNull(i) = max(cohN(Ub .* exp(-1j*PhiExc*b/NP)), cohN(Ub .* exp(1j*PhiExc*b/NP)));
end
H.Detect.Threshold = prctile(H.Detect.CohNull, 95);
H.Detect.Detected = H.Detect.CohExc > H.Detect.Threshold;
H.T_SBM = Info.Ts_SBM;
H.Time_equiv_ms = (0:nF-1) * Info.Ts_SBM * 1e3;
end

function c = nbCoh(U, dim)
% Normalised lag-1 neighbour coherence of a complex map along dim (NaN-safe).
if dim == 1, P = U(1:end-1, :) .* conj(U(2:end, :));
else,        P = U(:, 1:end-1) .* conj(U(:, 2:end)); end
c = abs(mean(P(:), 'omitnan')) / mean(abs(U(:)).^2, 'omitnan');
end
