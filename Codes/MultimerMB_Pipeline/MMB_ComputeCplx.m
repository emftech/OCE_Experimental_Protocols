function C = MMB_ComputeCplx(m, posIdx, DC_spectrum, hannK, tIdx, zIdx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_ComputeCplx
%
% Complex OCT signal for a subset of lateral positions, with exactly the
% same math as ReadCplx_Volume_BS (global DC spectrum subtraction, Hann
% window along the spectrum, FFT, depth crop, time crop), but reading the
% raw fringes from a memory-mapped file so that only one meridian is in
% memory at a time (a 32-meridian file would otherwise need ~33 GB).
%
% Inputs:
%   m           : memmapfile with field raw [samples x A-lines x positions]
%   posIdx      : lateral positions to read
%   DC_spectrum : [samples x 1] median spectrum
%   hannK       : [samples x 1] Hann window
%   tIdx, zIdx  : time (A-line) and depth indices to keep
% Output:
%   C : complex [numel(posIdx) x numel(zIdx) x numel(tIdx)]
%       (same layout as Cplx_Matrix in PP_MBmode.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = complex(zeros(numel(posIdx), numel(zIdx), numel(tIdx)));
for i = 1:numel(posIdx)
    fr = double(m.Data.raw(:, tIdx, posIdx(i))) - DC_spectrum;
    F = fft(hannK .* fr);
    C(i,:,:) = F(zIdx,:);
end

end
