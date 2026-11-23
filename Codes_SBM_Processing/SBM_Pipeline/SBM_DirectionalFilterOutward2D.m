function [UDF, Info] = SBM_DirectionalFilterOutward2D(U, x_mm, y_mm, Source_mm, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_DirectionalFilterOutward2D
%
% Outward (propagation-direction) filter for the full-field SBM harmonic
% phasor, the 2D counterpart of MMB_DirectionalFilterOutward. In MB mode
% every meridian passes through the needle and each half is filtered in
% the k-f plane. With a full raster the wave of interest travels outward
% from the needle in every direction, while reflections from the phantom
% boundaries travel inward.
%
% With the time convention of the phasor (U ~ exp(-j k.r) for a wave
% travelling along k) and the MATLAB 2D FFT kernel exp(-j 2 pi kappa.r),
% a component travelling in direction n has its energy at kappa.n < 0.
% The plane is split into NSectors angular sectors around the needle; for
% sector s (outward direction n_s) the half-plane mask
%     H_s(kappa) = 0.5 erfc( kappa.n_s / (sqrt(2) sigma) ),
%     sigma = SigmaKappaRel * kappa0      (smoothed, as in the MB filter)
% keeps the components travelling away from the needle. The filtered
% sector fields are blended with cos^2 angular weights that sum to 1:
%     U_out(r) = sum_s w_s(theta(r)) IFFT{ FFT{U} H_s }
%
% Inputs:
%   U          : [Ny x Nx] complex phasor, NaN outside the support
%   x_mm, y_mm : axes (y ascending)
%   Source_mm  : [x0 y0] needle tip
% Optional name-value:
%   'NSectors' (8), 'SigmaKappaRel' (0.5), 'PadFactor' (4),
%   'kappa0_cyc_per_mm' ([] = dominant ring of |FFT{U}|)
% Outputs:
%   UDF  : filtered phasor (NaN where U is NaN)
%   Info : kappa0_cyc_per_mm, OutwardEnergyFraction, Spec_Before/After
%          (|FFT| on kappa_x, kappa_y axes), Params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'NSectors', 8);
addParameter(p, 'SigmaKappaRel', 0.5);
addParameter(p, 'PadFactor', 4);
addParameter(p, 'kappa0_cyc_per_mm', []);
parse(p, varargin{:});
opts = p.Results;

[Ny, Nx] = size(U);
sup = isfinite(U);
U0 = U;  U0(~sup) = 0;
dx = mean(diff(x_mm));  dy = mean(diff(y_mm));
Px = 2^nextpow2(opts.PadFactor*Nx);  Py = 2^nextpow2(opts.PadFactor*Ny);
S = fft2(U0, Py, Px);
kx = (-floor(Px/2):ceil(Px/2)-1)/(Px*dx);       % cycles/mm (fftshifted)
ky = (-floor(Py/2):ceil(Py/2)-1)/(Py*dy);
[KX, KY] = meshgrid(kx, ky);
Ss = fftshift(S);

kap0 = opts.kappa0_cyc_per_mm;
if isempty(kap0)                                % dominant radius of |FFT|
    Kr = hypot(KX, KY);
    edges = linspace(0, max(Kr(:))/2, 200);
    [~, ~, bin] = histcounts(Kr(:), edges);
    A = accumarray(max(bin, 1), abs(Ss(:)).^2, [numel(edges)-1 1]);
    A(1:2) = 0;                                 % ignore DC
    [~, i] = max(A);
    kap0 = 0.5*(edges(i) + edges(i+1));
end
sig = opts.SigmaKappaRel*kap0;

[X, Y] = meshgrid(x_mm, y_mm);
th = atan2(Y - Source_mm(2), X - Source_mm(1));
N = opts.NSectors;  Dth = 2*pi/N;
UDF = zeros(Ny, Nx);
for s = 0:N-1
    ths = -pi + s*Dth;
    n = [cos(ths) sin(ths)];
    H = 0.5*erfc((KX*n(1) + KY*n(2))/(sqrt(2)*sig));
    Us = ifft2(ifftshift(Ss.*H));
    Us = Us(1:Ny, 1:Nx);
    d = angle(exp(1j*(th - ths)));              % wrapped angular distance
    w = cos(pi*d/(2*Dth)).^2 .* (abs(d) < Dth);
    UDF = UDF + w.*Us;
end
UDF(~sup) = NaN;

Info = struct();
Info.kappa0_cyc_per_mm = kap0;
Info.OutwardEnergyFraction = sum(abs(UDF(sup)).^2)/sum(abs(U0(sup)).^2);
Info.kx_cyc_per_mm = kx;  Info.ky_cyc_per_mm = ky;
Info.Spec_Before = single(abs(Ss));
SA = fft2(nan2zero(UDF), Py, Px);
Info.Spec_After = single(abs(fftshift(SA)));
Info.Params = opts;
end

function A = nan2zero(A)
A(~isfinite(A)) = 0;
end
