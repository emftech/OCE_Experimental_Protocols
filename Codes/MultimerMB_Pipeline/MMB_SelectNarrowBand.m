function Band = MMB_SelectNarrowBand(V, fs, f_exc, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_SelectNarrowBand
%
% Selects a NARROW temporal pass band around the excitation frequency from
% the measured motion spectrum of a Multimeridian MB-mode acquisition.
%
% Rule (07-10-2026 processing):
%   1) Mean power spectrum P(f) of all valid surface traces (mean removed,
%      Hann window, zero-padded FFT).
%   2) Search the spectral peak closest to f_exc inside
%      f_exc*(1 -/+ SearchRel). The peak is "clearly defined" when it is a
%      local maximum with prominence >= MinProminence_dB and both -3 dB
%      (half-power) points exist around it.
%   3) Pass-band width = -3 dB width of that peak, but never wider than
%      MaxBW_Hz (the -3 dB width is used "at the most"). If the peak is not
%      clearly defined, the band is centred at f_exc with width MaxBW_Hz.
%   4) The analysis frequency f_an (used for the phasor AND for c = 2 pi f/|k|)
%      is the peak frequency when the peak is clearly defined, otherwise
%      f_exc.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   V     : [nTraces x Nt] motion traces (rows with NaN are ignored)
%   fs    : sampling frequency [Hz]
%   f_exc : nominal excitation frequency [Hz]
%
% Optional name-value parameters:
%   'SearchRel'        : relative half-width of the peak search window (0.20)
%   'MaxBW_Hz'         : maximum pass-band width [Hz] (600)
%   'MinProminence_dB' : prominence for a clearly defined peak (3)
%   'Nfft'             : FFT length (2^15)
%
% Output struct Band:
%   f_exc, f_peak, f_an, PeakIsClear, f3dB_lo, f3dB_hi, BW3dB_Hz,
%   BW_Hz, f_lo, f_hi (final pass band), Prominence_dB,
%   freqs, P (mean power spectrum, linear), P_dB (normalized to the peak),
%   Resolution_Hz (= fs/Nt, native resolution of the record)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'SearchRel', 0.20);
addParameter(p, 'MaxBW_Hz', 600);
addParameter(p, 'MinProminence_dB', 3);
addParameter(p, 'Nfft', 2^15);
parse(p, varargin{:});
opts = p.Results;

%% ===================== Mean power spectrum ===============================

V = V(all(isfinite(V), 2), :);
Nt = size(V, 2);
V = (V - mean(V, 2)) .* hann(Nt).';

S = mean(abs(fft(V, opts.Nfft, 2)).^2, 1);
nHalf = floor(opts.Nfft/2);
freqs = (0:nHalf-1) * fs/opts.Nfft;
P = S(1:nHalf);

%% ===================== Peak nearest to f_exc =============================

win = freqs >= f_exc*(1-opts.SearchRel) & freqs <= f_exc*(1+opts.SearchRel);
PdB = 10*log10(P/max(P(win)));

[pk, locs, ~, prom] = findpeaks(PdB);
inWin = win(locs);
pk = pk(inWin); locs = locs(inWin); prom = prom(inWin);

Band = struct();
Band.f_exc = f_exc;
Band.PeakIsClear = false;
Band.f_peak = NaN;
Band.Prominence_dB = NaN;
Band.f3dB_lo = NaN;
Band.f3dB_hi = NaN;
Band.BW3dB_Hz = NaN;

if ~isempty(locs)
    % Highest peak in the window (ties -> closest to f_exc)
    [~, order] = sortrows([-pk(:), abs(freqs(locs(:)).' - f_exc)]);
    i0 = locs(order(1));
    Band.f_peak = freqs(i0);
    Band.Prominence_dB = prom(order(1));

    % -3 dB (half power) points around the peak
    lvl = PdB(i0) - 3;
    iL = find(PdB(1:i0) < lvl, 1, 'last');
    iR = i0 - 1 + find(PdB(i0:end) < lvl, 1, 'first');

    if ~isempty(iL) && ~isempty(iR)
        Band.f3dB_lo = interp1(PdB([iL iL+1]), freqs([iL iL+1]), lvl);
        Band.f3dB_hi = interp1(PdB([iR-1 iR]), freqs([iR-1 iR]), lvl);
        Band.BW3dB_Hz = Band.f3dB_hi - Band.f3dB_lo;
        Band.PeakIsClear = Band.Prominence_dB >= opts.MinProminence_dB;
    end
end

%% ===================== Final pass band ===================================

if Band.PeakIsClear
    Band.f_an = Band.f_peak;
    Band.BW_Hz = min(Band.BW3dB_Hz, opts.MaxBW_Hz);
else
    Band.f_an = f_exc;
    Band.BW_Hz = opts.MaxBW_Hz;
end

Band.f_lo = Band.f_an - Band.BW_Hz/2;
Band.f_hi = Band.f_an + Band.BW_Hz/2;

Band.freqs = freqs;
Band.P = P;
Band.P_dB = 10*log10(P/max(P(win)));
Band.Resolution_Hz = fs/Nt;
Band.Options = opts;

end
