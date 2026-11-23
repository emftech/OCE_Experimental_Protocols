function [Recovered, Phasor, PhiAcq, Raw] = SBM_TwinOCT_Recover(Frames, Info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_Recover
%
% Time-domain acquisition-phase compensation (thesis Sec. 5.2, Eq. 57-59)
% using the EXACT acquisition instant of every pixel, taken from the
% hardware timeline (Info), not from a nominal raster model.
%
% For pixel (b,a) in frame m the A-line is recorded at the global sample
%   n = s_1 + (m-1) F + tau(b,a),   F = K P + N_phi
% For a P-periodic (steady-state harmonic) field this is the equivalent
% time (m-1) N_phi Ts + (s_1 + tau(b,a)) Ts, so the frame series is a
% uniformly sampled waveform over N_P periods (Eq. 56) delayed by
%   phi_acq(b,a) = omega_exc (s_1 + tau(b,a)) Ts            (Eq. 57)
% With N_F*N_phi = N_P*P, f_exc falls exactly on DFT bin N_P, so
%   U_sampled = (2/N_F) sum_m u_m e^{-j 2 pi N_P (m-1)/N_F}
%   U_rec     = U_sampled e^{-j phi_acq}                     (Eq. 59)
% phi_acq is referenced to the excitation-sine phase (sin(2 pi n/P), n = 0
% at the first AO sample), so angle(U_rec) is the wave phase relative to
% the drive.
%
% Inputs:
%   Frames : [N_B x N_A x N_F] real displacement / phase (from
%            SBM_TwinOCT_Reorganize). Remove static phase first.
%   Info   : from SBM_TwinOCT_Waveform.
% Outputs:
%   Recovered : [N_B x N_A x N_F] real harmonic frames after compensation,
%               evaluated at the equivalent times (m-1) Ts_SBM.
%   Phasor    : complex phasor at f_exc after compensation (Eq. 59)
%   PhiAcq    : accumulated acquisition phase matrix [rad]
%   Raw       : complex phasor before compensation (Eq. 58)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nF = size(Frames,3);
NP = Info.Cfg.nPeriods;
m  = reshape(0:nF-1, 1, 1, []);
Raw = (2/nF) * sum(double(Frames) .* exp(-1j*2*pi*NP*m/nF), 3, 'omitnan');
PhiAcq = Info.PhiAcq;
Phasor = Raw .* exp(-1j*PhiAcq);
Recovered = real(Phasor .* exp(1j*2*pi*NP*m/nF));
end
