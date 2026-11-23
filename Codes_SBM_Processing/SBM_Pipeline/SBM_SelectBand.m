function Band = SBM_SelectBand(V, nTraces, f_exc, NP, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_SelectBand
%
% SBM counterpart of MMB_SelectNarrowBand: selects the temporal pass band
% from the mean spectrum of the SBM frame series. The frames are samples in
% EQUIVALENT time T_SBM, so bin k sits at k f_exc/N_P (f_exc = bin N_P) and
% the native resolution is f_exc/N_P (877 Hz at 285 us, N_P = 4).
%
% Rule (as in the MB processing of 07-10-2026):
%   1) mean power spectrum of all pixels (mean removed, Hann window,
%      zero-padded to Nfft);
%   2) highest local peak inside f_exc*(1 -/+ SearchRel) (ties -> closest
%      to f_exc). "Clearly defined" = prominence >= MinProminence_dB and
%      both -3 dB points exist;
%   3) pass-band width = its -3 dB width, at most MaxBW_Hz (default one bin,
%      f_exc/N_P: narrower than the native resolution is meaningless here).
%      If the peak is not clearly defined: centred at f_exc, width MaxBW_Hz;
%   4) f_an = peak frequency (UsePeakFrequency = true) or f_exc.
% Mode 'exact' skips all of it: f_an = f_exc, no filter (bin N_P only).
%
% CAUTION: only components locked to the excitation period (f_exc and its
% harmonics) sit at their true frequency in equivalent time. Anything else
% (e.g. the transducer resonance, table vibration) is aliased by the frame
% rate 1/(F Ts) (146 Hz at 285 us / 21 B-scans) onto an arbitrary bin, so
% a peak other than bin N_P must be checked with the phi_acq detection test.
%
%   Band = SBM_SelectBand(V, nTraces, f_exc, NP, opts)
% V       : frame series of the valid pixels, any shape with N_F along the
%           last dimension and nTraces elements per frame
% opts    : Mode ('exact' | 'peak'), SearchRel (0.5), MaxBW_Hz ([] = one
%           bin), MinProminence_dB (3), UsePeakFrequency (true), Nfft (1024)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def = struct('Mode', 'exact', 'SearchRel', 0.5, 'MaxBW_Hz', [], 'MinProminence_dB', 3, ...
    'UsePeakFrequency', true, 'Nfft', 1024);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end, end
df = f_exc / NP;
if isempty(opts.MaxBW_Hz), opts.MaxBW_Hz = df; end

Band = struct('Mode', opts.Mode, 'f_exc', f_exc, 'f_peak', NaN, 'PeakIsClear', false, ...
    'Prominence_dB', NaN, 'f3dB_lo', NaN, 'f3dB_hi', NaN, 'BW3dB_Hz', NaN, ...
    'f_an', f_exc, 'BW_Hz', NaN, 'Filtered', false, 'Resolution_Hz', df, 'Options', opts);
if strcmpi(opts.Mode, 'exact'), return; end

V = reshape(V, nTraces, []);
nF = size(V, 2);
fsEq = nF * df;
V = (V - mean(V, 2)) .* hann(nF).';
S = mean(abs(fft(V, opts.Nfft, 2)).^2, 1);
nHalf = floor(opts.Nfft/2);
freqs = (0:nHalf-1) * fsEq/opts.Nfft;
P = S(1:nHalf);

win = freqs >= f_exc*(1-opts.SearchRel) & freqs <= f_exc*(1+opts.SearchRel);
PdB = 10*log10(P/max(P(win)));
[pk, locs, ~, prom] = findpeaks(PdB);
inWin = win(locs);
pk = pk(inWin);  locs = locs(inWin);  prom = prom(inWin);
if ~isempty(locs)
    [~, order] = sortrows([-pk(:), abs(freqs(locs(:)).' - f_exc)]);
    i0 = locs(order(1));
    Band.f_peak = freqs(i0);
    if i0 > 1 && i0 < nHalf                         % parabolic refinement (zero-padding grid)
        y = PdB(i0-1:i0+1);  d = 0.5*(y(1) - y(3))/(y(1) - 2*y(2) + y(3));
        if isfinite(d) && abs(d) < 1, Band.f_peak = freqs(i0) + d*(freqs(2) - freqs(1)); end
    end
    Band.Prominence_dB = prom(order(1));
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
if Band.PeakIsClear
    Band.BW_Hz = min(Band.BW3dB_Hz, opts.MaxBW_Hz);
    if opts.UsePeakFrequency, Band.f_an = Band.f_peak; end
else
    Band.BW_Hz = opts.MaxBW_Hz;
end
Band.Filtered = true;
Band.f_lo = Band.f_an - Band.BW_Hz/2;
Band.f_hi = Band.f_an + Band.BW_Hz/2;
Band.freqs = freqs;
Band.P = P;
Band.P_dB = 10*log10(P/max(P));
end
