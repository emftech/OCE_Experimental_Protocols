function Meridian_Frames = OCE_MB_ReshapeToMeridians(SpaceTime, nPosPerMer, nMer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_MB_ReshapeToMeridians
%
% Converts the concatenated space-time map of a Multimeridian MB-mode
% acquisition into per-meridian space-time frames.
%
% Input:
%   SpaceTime : [nMer*nPosPerMer x Nt] matrix, e.g. loaded_phases_Border
%               from MotionAnalysisBorder(). Rows (b-1)*nPosPerMer+1 :
%               b*nPosPerMer belong to meridian b, in acquisition order.
%
% Output:
%   Meridian_Frames : [nPosPerMer x Nt x nMer], same layout used by the
%                     simulator (Meridian_Frames_MB).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nRows, Nt] = size(SpaceTime);

if nRows ~= nPosPerMer*nMer
    error(['SpaceTime has %d rows, but nPosPerMer*nMer = %d*%d = %d. ' ...
           'Check that SpaceTime is [lateral positions x time].'], ...
           nRows, nPosPerMer, nMer, nPosPerMer*nMer);
end

Meridian_Frames = permute(reshape(SpaceTime, nPosPerMer, nMer, Nt), [1 3 2]);

end
