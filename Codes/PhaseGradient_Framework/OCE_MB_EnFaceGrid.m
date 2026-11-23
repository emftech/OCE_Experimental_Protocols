function [x_mm, y_mm, Xq_mm, Yq_mm] = OCE_MB_EnFaceGrid(Geom, GridSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_MB_EnFaceGrid
%
% Regular en face grid covering the multimeridian field of view
% [-L/2, L/2] x [-L/2, L/2] around the pattern center.
%
% Uses MESHGRID convention: Xq_mm/Yq_mm are [Ny x Nx], rows = y.
%
% Inputs:
%   Geom     : struct from OCE_MB_MeridianGeometry
%   GridSize : scalar N (N x N grid) or [Ny Nx]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isscalar(GridSize)
    GridSize = [GridSize GridSize];
end

halfFOV = Geom.scanLength_mm/2;

x_mm = Geom.Center_mm(1) + linspace(-halfFOV, halfFOV, GridSize(2));
y_mm = Geom.Center_mm(2) + linspace(-halfFOV, halfFOV, GridSize(1));

[Xq_mm, Yq_mm] = meshgrid(x_mm, y_mm);

end
