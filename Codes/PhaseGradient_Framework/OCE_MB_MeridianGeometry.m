function Geom = OCE_MB_MeridianGeometry(nPosPerMer, nMer, scanLength_mm, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_MB_MeridianGeometry
%
% Rebuilds the en face (x,y) position of every lateral position acquired in
% a real Multimeridian MB-mode acquisition, in ACQUISITION ORDER.
%
% It follows the same convention as AcqVectors_MultimeridianMB.m (simulator)
% and IDXSelection.m (experimental pipeline):
%
%   - nMer meridians with angles theta = 0 : 180/nMer : 180*(1 - 1/nMer)
%   - each meridian spans the full scan length, centered on the excitation
%     point: r = linspace(-L/2, +L/2, nPosPerMer)
%   - every second meridian is acquired in reverse radial order
%
%       x = xc + r cos(theta),   y = yc + r sin(theta)
%
% Row i of meridian b therefore corresponds to row i of block b in the
% concatenated space-time map produced by MotionAnalysisBorder().
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   nPosPerMer    : lateral positions per meridian (ScanInfo.Bframes_in_3Dscan)
%   nMer          : number of meridians (ScanInfo.No_3Dscans)
%   scanLength_mm : meridian length [mm] (ScanInfo.Ver_scan_length_mm)
%
% Optional name-value parameters:
%   'ReverseAlternate' : every 2nd meridian acquired in reverse (default true)
%   'AngleOffset_deg'  : rotation added to all meridian angles (default 0)
%   'FlipRadial'       : invert the radial direction of ALL meridians,
%                        i.e. first sample at +L/2 instead of -L/2
%                        (default false)
%   'Center_mm'        : [xc yc] excitation / pattern center (default [0 0])
%
% Output struct Geom:
%   X_mm, Y_mm  : [nPosPerMer x nMer] en face position of each sample
%   s_mm        : [nPosPerMer x nMer] signed coordinate along the meridian
%                 (monotonic in acquisition order, used by the 1D PG fit)
%   r_mm        : [nPosPerMer x nMer] radial distance to the center
%   theta_rad   : [1 x nMer] meridian angles
%   ux, uy      : [1 x nMer] unit vector of increasing s for each meridian
%   scanLength_mm, nPosPerMer, nMer, Center_mm, ds_mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.FunctionName = 'OCE_MB_MeridianGeometry';
addParameter(p, 'ReverseAlternate', true);
addParameter(p, 'AngleOffset_deg', 0);
addParameter(p, 'FlipRadial', false);
addParameter(p, 'Center_mm', [0 0]);
parse(p, varargin{:});
opts = p.Results;

%% ===================== Radial sampling and angles ========================

r = linspace(-scanLength_mm/2, scanLength_mm/2, nPosPerMer).';

if opts.FlipRadial
    r = flipud(r);
end

theta = linspace(0, pi, nMer + 1);
theta = theta(1:end-1) + deg2rad(opts.AngleOffset_deg);

rr = repmat(r, 1, nMer);

% Match multimeridian acquisition convention:
% every second meridian is acquired in reverse radial order.
if opts.ReverseAlternate
    rr(:,2:2:end) = rr(end:-1:1,2:2:end);
end

tt = repmat(theta, nPosPerMer, 1);

%% ===================== En face coordinates ===============================

Geom = struct();

Geom.X_mm = opts.Center_mm(1) + rr .* cos(tt);
Geom.Y_mm = opts.Center_mm(2) + rr .* sin(tt);

% Signed coordinate along the meridian. It is kept in acquisition order,
% so reversed meridians have decreasing s. The 1D PG fit uses s directly,
% therefore the sign of k along reversed meridians is flipped, which does
% not affect |k|.
Geom.s_mm = rr;
Geom.r_mm = abs(rr);

Geom.theta_rad = theta;
Geom.angles_deg = rad2deg(theta);

Geom.ux = cos(theta);
Geom.uy = sin(theta);

Geom.scanLength_mm = scanLength_mm;
Geom.nPosPerMer = nPosPerMer;
Geom.nMer = nMer;
Geom.Center_mm = opts.Center_mm;
Geom.ds_mm = scanLength_mm / (nPosPerMer - 1);
Geom.Options = opts;

end
