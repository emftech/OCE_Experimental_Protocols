function cmap = MMB_RedBlue(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 20-10-2026
%
% MMB_RedBlue
%
% Diverging blue - white - red colormap with white exactly at zero (odd n),
% for signed wave fields (particle velocity). Negative = blue
% [0.13 0.40 0.75], zero = white, positive = red [0.70 0.09 0.17].
% In the depth videos white is never drawn: the overlay opacity goes to 0
% at v = 0, so the B-mode shows through (see MMB_WriteDepthVideo).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1, n = 255; end
if mod(n, 2) == 0, n = n + 1; end

t = linspace(0, 1, (n+1)/2)';
blue = [0.13 0.40 0.75];
red  = [0.70 0.09 0.17];
lower = blue + t.*(1 - blue);                 % blue -> white
upper = 1 + t(2:end).*(red - 1);              % white -> red
cmap = [lower; upper];

end
