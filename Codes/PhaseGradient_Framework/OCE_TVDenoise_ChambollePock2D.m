function [U, Info] = OCE_TVDenoise_ChambollePock2D(F, mu, Mask, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_TVDenoise_ChambollePock2D
%
% Solves the PG-TV inverse problem (thesis Eqs. 17-18) on the pixels of
% Mask:
%
%       U = argmin_U 0.5 ||U - F||_2^2 + mu * TV(U)
%
%       TV(U) = sum_ij sqrt( (dU/dx)_ij^2 + (dU/dy)_ij^2 )   (isotropic)
%
% using the Chambolle-Pock primal-dual algorithm.
%
% Boundary handling ('Boundary' option):
%
%   'neumann' (default):
%       Finite differences are only taken between two pixels that are BOTH
%       inside Mask. Pixels outside Mask are fully decoupled, so the TV
%       term never "sees" an artificial jump at the edge of the valid
%       support (natural Neumann condition on the mask boundary).
%
%   'legacy':
%       Reproduces TVDenoise_ChambollePock2D.m from the simulation study:
%       pixels outside Mask are fixed to 0 and the gradient is taken over
%       the full grid. Use it only to reproduce thesis results; it biases
%       the regularized map towards 0 next to the mask boundary.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   F    : 2D map to regularize (e.g. raw wave-number map |k| [rad/m])
%   mu   : TV regularization parameter (same units as F)
%   Mask : logical support where the problem is solved
%
% Optional name-value parameters:
%   'NIterations' : default 500
%   'Tau'         : primal step, default 0.25
%   'Sigma'       : dual step,   default 0.25   (Tau*Sigma*8 <= 1)
%   'Theta'       : over-relaxation, default 1.0
%   'Boundary'    : 'neumann' (default) or 'legacy'
%
% Outputs:
%   U    : regularized map, NaN outside Mask
%   Info : struct with the options used and the final relative update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.FunctionName = 'OCE_TVDenoise_ChambollePock2D';
addParameter(p, 'NIterations', 500);
addParameter(p, 'Tau', 0.25);
addParameter(p, 'Sigma', 0.25);
addParameter(p, 'Theta', 1.0);
addParameter(p, 'Boundary', 'neumann', @(s) any(strcmpi(s, {'neumann','legacy'})));
parse(p, varargin{:});
opts = p.Results;

Mask = logical(Mask) & isfinite(F);

Info = opts;
Info.mu = mu;
Info.ValidPixels = nnz(Mask);

if mu == 0 || ~any(Mask(:))
    U = F;
    U(~Mask) = NaN;
    Info.RelativeUpdate = 0;
    return;
end

if opts.Tau*opts.Sigma*8 > 1
    warning('Tau*Sigma*8 = %.3g > 1: Chambolle-Pock may not converge.', ...
        opts.Tau*opts.Sigma*8);
end

%% ===================== Edge masks ========================================

[Nx, Ny] = size(F);

switch lower(opts.Boundary)

    case 'neumann'
        % Only differences between two valid pixels are penalized.
        Mx = false(Nx, Ny);
        My = false(Nx, Ny);
        Mx(1:end-1,:) = Mask(1:end-1,:) & Mask(2:end,:);
        My(:,1:end-1) = Mask(:,1:end-1) & Mask(:,2:end);

        F0 = F;
        F0(~Mask) = 0;

    case 'legacy'
        Mx = true(Nx, Ny);
        My = true(Nx, Ny);
        Mx(end,:) = false;
        My(:,end) = false;

        F0 = F;
        F0(~Mask) = 0;
end

%% ===================== Chambolle-Pock iterations =========================

U = F0;
Ubar = U;

Px = zeros(Nx, Ny);
Py = zeros(Nx, Ny);

tau = opts.Tau;
sigma = opts.Sigma;
theta = opts.Theta;

for it = 1:opts.NIterations

    % Dual ascent + projection onto the mu-ball
    [Gx, Gy] = maskedGradient(Ubar, Mx, My);

    Px = Px + sigma*Gx;
    Py = Py + sigma*Gy;

    NormP = max(1, sqrt(Px.^2 + Py.^2)/mu);

    Px = Px ./ NormP;
    Py = Py ./ NormP;

    % Primal descent (proximal operator of the quadratic data term)
    Uold = U;

    U = (U + tau*divergence2D(Px, Py) + tau*F0) / (1 + tau);

    U(~Mask) = F0(~Mask);

    Ubar = U + theta*(U - Uold);
end

Info.RelativeUpdate = norm(U(Mask) - Uold(Mask)) / max(norm(U(Mask)), eps);

U(~Mask) = NaN;

end

%% ======================================================================= %
%                          LOCAL FUNCTIONS
% ======================================================================== %

function [Gx, Gy] = maskedGradient(U, Mx, My)

Gx = zeros(size(U));
Gy = zeros(size(U));

Gx(1:end-1,:) = U(2:end,:) - U(1:end-1,:);
Gy(:,1:end-1) = U(:,2:end) - U(:,1:end-1);

Gx(~Mx) = 0;
Gy(~My) = 0;

end

function Div = divergence2D(Px, Py)
% Negative adjoint of the forward-difference gradient. Px and Py are zero
% on edges that are not part of the masked graph, so this is also the
% adjoint of the masked gradient.

Div = zeros(size(Px));

Div(1,:) = Px(1,:);
Div(2:end-1,:) = Px(2:end-1,:) - Px(1:end-2,:);
Div(end,:) = -Px(end-1,:);

Div(:,1) = Div(:,1) + Py(:,1);
Div(:,2:end-1) = Div(:,2:end-1) + Py(:,2:end-1) - Py(:,1:end-2);
Div(:,end) = Div(:,end) - Py(:,end-1);

end
