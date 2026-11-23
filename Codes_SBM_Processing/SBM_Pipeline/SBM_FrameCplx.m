function C = SBM_FrameCplx(m, f, DC_spectrum, hannK, zIdx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_FrameCplx
%
% Complex OCT signal of ONE SBM frame (all C*R A-lines of the frame, in
% the order the ATS9360 stored them), with exactly the math of
% MMB_ComputeCplx / ReadCplx_Volume_BS: global DC spectrum subtraction,
% Hann window along the spectrum, FFT, depth crop.
%
% Inputs:
%   m           : memmapfile with field raw [samples x (R*C) x N_F]
%   f           : frame index
%   DC_spectrum : [samples x 1] median spectrum
%   hannK       : [samples x 1] Hann window
%   zIdx        : depth indices to keep
% Output:
%   C : complex [numel(zIdx) x (R*C)] (depth x A-line in record order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = fft(hannK .* (double(m.Data.raw(:,:,f)) - DC_spectrum));
C = F(zIdx, :);
end
