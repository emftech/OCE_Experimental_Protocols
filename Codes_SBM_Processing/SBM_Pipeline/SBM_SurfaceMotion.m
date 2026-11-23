function V = SBM_SurfaceMotion(Band, Valid, Options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_SurfaceMotion
%
% Surface motion of every raster position across the N_F SBM frames, the
% SBM counterpart of MotionAnalysisBorder.m (PP_MBmode). In MB mode the
% time axis is the A-lines of one M-scan; in SBM it is the FRAME index,
% i.e. the equivalent-time samples spaced T_SBM = N_phi*Ts (thesis Eq. 44).
%
% Options.MotionType = 2 (particle velocity, default, as in MB):
%   Loupas 2D estimator (Loupas2D_Fast, axial window LoupasAxialWin) at
%   lag 1 between consecutive frames, averaged over Delta+1 depth rows
%   starting LoupasAxialWin/2 above the surface (exactly the depth
%   placement of MotionAnalysisBorder). The N_F frames span exactly N_P
%   excitation periods, so the lag is taken CIRCULARLY (frame N_F -> 1):
%   the output keeps N_F samples and the excitation stays on DFT bin N_P.
%   Output [rad per equivalent step]; multiply by lambda0/(4 pi n)/T_SBM
%   for mm/s.
% Options.MotionType = 1 (displacement):
%   -angle( sum_depth s_m conj(<s>_m) ) over rows Border..Border+Delta,
%   i.e. the phase relative to the frame-averaged (static) signal [rad].
%   Multiply by lambda0/(4 pi n) for displacement.
% Options.Estimator = 'band' (MotionType 2): phase of the complex lag-1
%   correlation summed over the BandDepth_px rows below the surface,
%   -angle(sum_z s_{m+1} conj(s_m)) (Loupas sign convention). On the real
%   SBM data it is ~4x less noisy than the 10-px Loupas window, half of
%   which lies above the surface (18-10-2026: 0.21 vs 0.9 rad after bulk
%   compensation); default 'loupas' keeps the MB estimator.
% Options.SmoothingWinPer : lowess smoothing along frames ([] = none,
%   recommended: the harmonic is extracted from one exact DFT bin).
%
% Inputs:
%   Band  : complex [nRows x Nper x N_F] depth band around the surface,
%           row 1 = Border - round(LoupasAxialWin/2) (SBM_SurfaceBand)
%   Valid : logical [Nper x 1] (surface found and band inside the image)
% Output:
%   V : [Nper x N_F] motion (NaN for invalid A-lines)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, Nper, nF] = size(Band);
M = Options.LoupasAxialWin;  Delta = Options.Delta;
h = round(M/2);
V = nan(Nper, nF);
useBand = isfield(Options, 'Estimator') && strcmpi(Options.Estimator, 'band');
if useBand && Options.MotionType == 2
    S = double(Band(:, Valid, :));
    V(Valid, :) = -reshape(angle(sum(S(:, :, [2:nF 1]) .* conj(S), 1)), [], nF);
    return
end
for j = find(Valid(:))'
    S = reshape(Band(:, j, :), [], nF);
    if Options.MotionType == 2
        D = Loupas2D_Fast(double(S(:, [1:nF 1])), M);      % circular lag 1
        v = mean(D(1:Delta+1, :), 1);
    else
        rows = h + (1:Delta+1);                            % Border .. Border+Delta
        Sb = double(S(rows, :));
        v = -angle(sum(Sb .* conj(mean(Sb, 2)), 1));
        v = unwrap(v);
    end
    if ~isempty(Options.SmoothingWinPer)
        v = smooth(v, Options.SmoothingWinPer, 'lowess')';
    end
    V(j, :) = v;
end
end
