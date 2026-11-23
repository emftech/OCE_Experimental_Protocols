function [Phasor, Info] = OCE_HarmonicPhasor(Data, dimT, fs, f_exc, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_HarmonicPhasor
%
% Isolates the harmonic component of the measured motion at the excitation
% frequency (thesis Eqs. 9-10):
%
%       V(x,y,w0) = F_t{ v(x,y,t) } |_(w = w0) = A(x,y) exp(j theta(x,y))
%
% Two evaluation methods are available:
%
%   'dtft'   : exact single-frequency DTFT at f_exc (default).
%              V = sum_n w[n] (v[n] - mean(v)) exp(-j 2 pi f_exc n / fs)
%              No bin-mismatch error, recommended for experimental data
%              where the record holds only a few excitation periods.
%
%   'fftbin' : zero-padded FFT and nearest bin to f_exc, as in the
%              simulation codes (PhaseGradient2D_FromSims_v3). Kept for
%              reproducibility of the thesis workflow.
%
% Traces that are partially NaN are filled by linear interpolation along
% time; traces that are entirely NaN give a NaN phasor.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   Data  : N-D real array with time along dimension dimT
%   dimT  : temporal dimension of Data
%   fs    : temporal sampling frequency [Hz] (A-line rate)
%   f_exc : excitation frequency [Hz]
%
% Optional name-value parameters:
%   'Method'       : 'dtft' (default) or 'fftbin'
%   'TemporalHann' : apply a Hann window before the transform (default true)
%   'RemoveMean'   : remove the temporal mean of each trace (default true)
%   'Nfft'         : FFT length for 'fftbin' (default 2^nextpow2(max(Nt,1024)))
%
% Outputs:
%   Phasor : complex array, same size as Data with dimT removed (squeezed
%            to the remaining dimensions in their original order)
%   Info   : struct with f_used, Method, Nt, fs, f_exc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.FunctionName = 'OCE_HarmonicPhasor';
addParameter(p, 'Method', 'dtft', @(s) any(strcmpi(s, {'dtft','fftbin'})));
addParameter(p, 'TemporalHann', true);
addParameter(p, 'RemoveMean', true);
addParameter(p, 'Nfft', []);
parse(p, varargin{:});
opts = p.Results;

%% ===================== Move time to the first dimension ==================

sz = size(Data);
nd = max(ndims(Data), dimT);
sz(end+1:nd) = 1;

perm = [dimT, setdiff(1:nd, dimT)];
X = permute(double(Data), perm);

Nt = sz(dimT);
spatialSize = sz(perm(2:end));
X = reshape(X, Nt, []);

%% ===================== Handle NaN samples ================================

allNaN = all(~isfinite(X), 1);
partNaN = any(~isfinite(X), 1) & ~allNaN;

tIdx = (1:Nt).';

for c = find(partNaN)
    trace = X(:,c);
    valid = isfinite(trace);
    if nnz(valid) == 1
        X(:,c) = trace(valid);
    else
        X(~valid,c) = interp1(tIdx(valid), trace(valid), tIdx(~valid), 'linear', 'extrap');
    end
end

X(:,allNaN) = 0;

%% ===================== Mean removal and temporal window ==================

if opts.RemoveMean
    X = X - mean(X, 1);
end

if opts.TemporalHann
    X = X .* hann(Nt);
end

%% ===================== Harmonic component at f_exc =======================

switch lower(opts.Method)

    case 'dtft'
        n = (0:Nt-1);
        kernel = exp(-1j*2*pi*f_exc*n/fs);
        P = kernel * X;
        f_used = f_exc;

    case 'fftbin'
        Nfft = opts.Nfft;
        if isempty(Nfft)
            Nfft = 2^nextpow2(max(Nt,1024));
        end
        freqs = (0:Nfft-1) * (fs/Nfft);
        [~, kpos] = min(abs(freqs - f_exc));
        n = (0:Nt-1);
        kernel = exp(-1j*2*pi*(kpos-1)*n/Nfft);
        P = kernel * X;
        f_used = freqs(kpos);
end

P(allNaN) = NaN;

if numel(spatialSize) == 1
    spatialSize = [spatialSize 1];
end

Phasor = reshape(P, spatialSize);

Info = struct();
Info.Method = lower(opts.Method);
Info.f_exc = f_exc;
Info.f_used = f_used;
Info.fs = fs;
Info.Nt = Nt;
Info.RecordLength_s = Nt/fs;
Info.NumExcitationPeriods = Nt*f_exc/fs;
Info.TemporalHann = opts.TemporalHann;

end
