function MMB_SaveFigure(fig, base, resolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_SaveFigure
%
% Saves a figure as .fig (editable in MATLAB) and .png (slides/notes).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    resolution = 200;
end

folder = fileparts(base);
if ~isempty(folder) && ~exist(folder, 'dir')
    mkdir(folder);
end

savefig(fig, [base '.fig']);
exportgraphics(fig, [base '.png'], 'Resolution', resolution);

end
