function [Band, Valid] = SBM_SurfaceBand(C, BorderIdx, Options, Band, f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_SurfaceBand
%
% Extracts, from the complex signal of one frame, the depth band that the
% Loupas surface estimator needs for every A-line:
%   rows Border - round(M/2) ... Border - round(M/2) + Delta + M - 1
% (M = LoupasAxialWin), or, for Options.Estimator = 'band', rows
%   Border ... Border + BandDepth_px - 1,
% so that only nRows x Nper x N_F complex values are
% kept in memory instead of the whole volume.
%
%   [Band, Valid] = SBM_SurfaceBand(C, BorderIdx, Options, Band, f)
% C         : [nZ x Nper] complex frame (SBM_FrameCplx)
% BorderIdx : [1 x Nper] surface depth index (in the cropped depth axis)
% Band      : preallocated single complex [nRows x Nper x N_F] ([] -> create)
% f         : frame index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = Options.LoupasAxialWin;  Delta = Options.Delta;
[nZ, Nper] = size(C);
if isfield(Options, 'Estimator') && strcmpi(Options.Estimator, 'band')
    nRows = Options.BandDepth_px;                    % rows Border .. Border+BandDepth-1
    r0 = round(BorderIdx(:)');
else
    nRows = Delta + M;
    r0 = round(BorderIdx(:)') - round(M/2);
end
Valid = isfinite(r0) & r0 >= 1 & r0 + nRows - 1 <= nZ;
if isempty(Band), Band = complex(zeros(nRows, Nper, 0, 'single')); end
rows = r0 + (0:nRows-1)';                       % [nRows x Nper]
rows(:, ~Valid) = 1;
idx = sub2ind([nZ Nper], rows, repmat(1:Nper, nRows, 1));
Band(:, :, f) = single(C(idx));
end
