function [bRep, Info] = MMB_SelectRepresentativeMeridian(VideoMask, Bmode_IntLog, Geom, Source_mm, MinBmode_dB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% MMB_SelectRepresentativeMeridian
%
% Chooses the meridian shown in the single-meridian depth video. The score
% of each meridian is the number of pixels where the depth velocity can be
% trusted (inside the tissue mask, needle columns removed, smoothed B-mode
% >= MinBmode_dB). The meridian with the highest score is chosen; among
% meridians within 5 % of the best score, the one passing closest to the
% needle tip is preferred (its space-time "V" is the most symmetric).
%
% Inputs:
%   VideoMask    : [nZ x nPos*nMer] logical (tissue, no needle columns)
%   Bmode_IntLog : [nZ x nPos*nMer] B-mode [dB]
%   Geom         : OCE_MB_MeridianGeometry struct
%   Source_mm    : [x0 y0] needle tip
%   MinBmode_dB  : reliability threshold (default 90)
% Outputs:
%   bRep : meridian index;  Info : Score, DistToTip_mm, Angle_deg, s0_mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5, MinBmode_dB = 90; end
nPos = Geom.nPosPerMer; nMer = Geom.nMer;
Rel = VideoMask & medfilt2(Bmode_IntLog, [5 5], 'symmetric') >= MinBmode_dB;
score = zeros(1, nMer);
for b = 1:nMer
    score(b) = nnz(Rel(:, (b-1)*nPos + (1:nPos)));
end
dTip = abs(Source_mm(1)*Geom.uy - Source_mm(2)*Geom.ux);     % distance tip - meridian line
cand = find(score >= 0.95*max(score));
[~, k] = min(dTip(cand));
bRep = cand(k);

Info = struct('Score', score, 'DistToTip_mm', dTip, 'Angle_deg', Geom.angles_deg, ...
    's0_mm', Source_mm(1)*Geom.ux(bRep) + Source_mm(2)*Geom.uy(bRep), 'Meridian', bRep);

end
