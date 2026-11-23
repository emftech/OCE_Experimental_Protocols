function [Vc, Bulk] = SBM_BulkCompensation(V, Window, Scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_BulkCompensation
%
% Removes the bulk (whole-sample) motion from the lag-1 SBM surface phase.
% Consecutive SBM frames are F*Ts apart (6.9 ms for 285 us / 21 B-scans),
% so any slow vibration of the table, phantom or optics (tens to hundreds
% of Hz, nm amplitude) changes the phase between frames by up to +-pi,
% almost equally for all the A-lines acquired close in time. That bulk term
% is aliased all over the frame spectrum and buries bin N_P (first real
% acquisitions 18-10-2026: 1.8 rad -> 0.16 rad after compensation).
%
% The bulk phase of A-line n (record order = acquisition order) is the
% circular moving mean of the lag-1 phase over Window A-lines centred on n
% (Window = P = one raster line = one excitation period, so a locked
% excitation is not averaged against itself in time), and it is subtracted.
%
%   [Vc, Bulk] = SBM_BulkCompensation(V, Window, Scale)
% V      : [Nper x N_F] lag-1 motion in record order (NaN = invalid A-line)
% Window : moving-window length in A-lines (e.g. Info.P)
% Scale  : units of V per radian (e.g. VelScale for mm/s; 1 for rad)
% Vc     : compensated motion, same units;  Bulk : removed bulk [rad]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ph = V / Scale;
E = exp(1j*Ph);
E(~isfinite(Ph)) = 0;
W = double(isfinite(Ph));
Bulk = angle(movsum(E, Window, 1) ./ max(movsum(W, Window, 1), 1));
Vc = angle(exp(1j*(Ph - Bulk))) * Scale;
Vc(~isfinite(V)) = NaN;
end
