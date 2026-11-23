function [f_dom, Info] = OCE_EstimateDominantFrequency(SpaceTime, fs, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_EstimateDominantFrequency
%
% Frequency of the maximum of the mean power spectrum of all valid traces
% (Hann window, heavy zero padding). Used as a sanity check of the
% excitation frequency: the phase speed c = 2*pi*f / |k| scales linearly
% with the frequency used by the phase-gradient estimator.
%
% Inputs:
%   SpaceTime : [nTraces x Nt] motion traces (NaN rows are ignored)
%   fs        : sampling frequency [Hz]
%
% Optional name-value parameters:
%   'Band' : [fmin fmax] search band in Hz (default [50 fs/4])
%   'Nfft' : default 2^16
%
% Outputs:
%   f_dom : dominant frequency [Hz]
%   Info  : struct with freqs, meanPower, resolution (fs/Nt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'Band', [50 fs/4]);
addParameter(p, 'Nfft', 2^16);
parse(p, varargin{:});
opts = p.Results;

V = SpaceTime(all(isfinite(SpaceTime), 2), :);
Nt = size(V, 2);

V = (V - mean(V, 2)) .* hann(Nt).';

P = mean(abs(fft(V, opts.Nfft, 2)).^2, 1);
freqs = (0:opts.Nfft-1) * fs/opts.Nfft;

inBand = freqs >= opts.Band(1) & freqs <= opts.Band(2);
Pb = P;
Pb(~inBand) = 0;

[~, imax] = max(Pb);
f_dom = freqs(imax);

Info = struct();
Info.freqs = freqs(1:floor(opts.Nfft/2));
Info.meanPower = P(1:floor(opts.Nfft/2));
Info.Resolution_Hz = fs/Nt;
Info.nTraces = size(V, 1);

end
