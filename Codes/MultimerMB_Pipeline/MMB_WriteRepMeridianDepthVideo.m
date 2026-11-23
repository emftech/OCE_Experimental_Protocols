function CL = MMB_WriteRepMeridianDepthVideo(VelMer, BmodeMer, MaskMer, Geom, b, Source_mm, Zaxis, Time_ms, outFile, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% MMB_WriteRepMeridianDepthVideo
%
% Depth motion video of ONE meridian (the representative one, see
% MMB_SelectRepresentativeMeridian), drawn along its geometric coordinate
% s (columns of reversed meridians are re-ordered so that s increases from
% left to right) with the projection of the needle tip s0 marked. Same
% overlay as MMB_WriteDepthVideo (blue-white-red, transparent at 0,
% limits = P90 of |v| over reliable pixels of this meridian).
%
% Inputs:
%   VelMer   : [nZ x nPos x Nt] depth velocity of meridian b (acquisition order)
%   BmodeMer, MaskMer : [nZ x nPos] B-mode [dB] and tissue mask of meridian b
%   Geom, b, Source_mm : geometry, meridian index, needle tip [x0 y0]
%   Zaxis, Time_ms, outFile : as in MMB_WriteDepthVideo
%   varargin : passed to MMB_WriteDepthVideo (FrameIdx, IThresh, ...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s, o] = sort(Geom.s_mm(:, b));
s0 = Source_mm(1)*Geom.ux(b) + Source_mm(2)*Geom.uy(b);

CL = MMB_WriteDepthVideo(VelMer(:, o, :), BmodeMer(:, o), MaskMer(:, o), s, Zaxis, Time_ms, outFile, ...
    'FigWidth', 1000, 'XLabel', sprintf('s (mm) along meridian %d (\\theta = %.1f°)', b, Geom.angles_deg(b)), ...
    'VLines', s0, 'VLineLabel', 'needle tip', varargin{:});

end
