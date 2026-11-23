function [Frames, x_mm, y_mm] = SBM_ToFrames(Data, Info, xW_mm, yW_mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_ToFrames
%
% Sorts per-A-line SBM data (record order, as stored by the ATS9360) into
% spatial frames on the raster grid, using the pixel map of the sidecar
% (Info.PixelIdx: (B-scan, column) -> A-line inside the frame record).
% The serpentine reversal is undone and the output follows the MESHGRID
% convention of the PG framework: rows = y ASCENDING (row 1 = bottom of the
% raster = Vertical Min), columns = x ascending (Horizontal Min -> Max).
%
% Inputs:
%   Data  : [Nper x N_F] (one value per A-line of the frame record and per
%           frame) or [Nz x Nper] (one record, e.g. B-mode)
%   Info  : sidecar Info (PixelIdx, Cfg.nAscans, Cfg.nBscans, R, C)
%   xW_mm, yW_mm : scan lengths (header / file name)
% Outputs:
%   Frames : [N_B x N_A x N_F]  or  [Nz x N_B x N_A]
%   x_mm   : [1 x N_A], y_mm : [1 x N_B] (ascending)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nA = Info.Cfg.nAscans;  nB = Info.Cfg.nBscans;  Nper = Info.R*Info.C;
PI = flipud(Info.PixelIdx);                    % row 1 = bottom (y ascending)
if size(Data, 1) == Nper                       % [Nper x N_F]
    nF = size(Data, 2);
    Frames = reshape(Data(PI(:), :), nB, nA, nF);
else                                           % [Nz x Nper]
    Frames = reshape(Data(:, PI(:)), size(Data, 1), nB, nA);
end
x_mm = linspace(-xW_mm/2, xW_mm/2, nA);
y_mm = linspace(-yW_mm/2, yW_mm/2, nB);
end
