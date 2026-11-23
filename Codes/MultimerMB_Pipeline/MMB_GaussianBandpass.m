function Y = MMB_GaussianBandpass(X, fs, fc, BW, dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_GaussianBandpass
%
% Zero-phase, narrowband Gaussian band-pass filter applied in the frequency
% domain along dimension dim:
%
%   H(f) = exp(-(|f| - fc)^2 / (2 sigma^2)),  sigma = (BW/2)/sqrt(ln 2)
%
% so that |H|^2 = 1/2 (-3 dB) exactly at fc -/+ BW/2. H is real and even,
% therefore the filter introduces no phase shift (same purpose as filtfilt).
% A Gaussian response has the shortest possible impulse response for a
% given bandwidth, which limits ringing at the edges of the short MB-mode
% records (449 samples = 2.2 ms at 200 kHz). Each trace is mean-removed and
% zero-padded to at least 3x its length before the FFT to avoid circular
% wrap-around.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   X   : real array; traces along dimension dim (default 2)
%   fs  : sampling frequency [Hz]
%   fc  : center frequency [Hz]
%   BW  : -3 dB bandwidth [Hz]
%   dim : time dimension (default 2)
% Output:
%   Y   : filtered array, same size as X. All-NaN traces stay NaN; isolated
%         NaN samples are treated as zeros.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5 || isempty(dim)
    dim = 2;
end

sz = size(X);
perm = [dim, setdiff(1:ndims(X), dim)];
Xp = permute(X, perm);
Nt = size(Xp, 1);
Xp = reshape(Xp, Nt, []);

allNaN = all(~isfinite(Xp), 1);
Xp(~isfinite(Xp)) = 0;
Xp = Xp - mean(Xp, 1);

Nfft = 2^nextpow2(3*Nt);
f = (0:Nfft-1).' * fs/Nfft;
f(f > fs/2) = f(f > fs/2) - fs;          % signed frequency axis

sigma = (BW/2)/sqrt(log(2));
H = exp(-(abs(f) - fc).^2 / (2*sigma^2));

Yp = real(ifft(fft(Xp, Nfft, 1) .* H, [], 1));
Yp = Yp(1:Nt, :);
Yp(:, allNaN) = NaN;

Yp = reshape(Yp, [Nt, sz(perm(2:end))]);
Y = ipermute(Yp, perm);

if isa(X, 'single')
    Y = single(Y);
end

end
