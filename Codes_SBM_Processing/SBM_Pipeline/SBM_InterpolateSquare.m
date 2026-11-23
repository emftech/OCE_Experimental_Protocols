function [Uq, xq, yq, Info] = SBM_InterpolateSquare(U, x_mm, y_mm, method, Nmin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_InterpolateSquare
%
% Experimental counterpart of Prepare_SBM_2DPhaseGradient (thesis SBM
% simulations): the recovered SBM field is interpolated onto a SQUARE grid
% of Nout x Nout pixels, Nout = max(N_A, N_B, Nmin) (Nmin = 100, the
% minimum grid used for the fair comparison with Multimeridian MB-mode),
% with the same interpolation method ('spline' in the thesis runs) and no
% smoothing. Works on a complex phasor [N_B x N_A] or a real movie
% [N_B x N_A x N_t] (frame by frame). NaN pixels (needle) are kept NaN:
% the field is filled before interpolation and a linearly interpolated
% validity mask removes every output pixel that touches them.
%
% Inputs:
%   U          : [N_B x N_A] complex or [N_B x N_A x N_t] real
%   x_mm, y_mm : axes (y ascending)
%   method     : 'spline' (default) | 'linear' | 'pchip' | 'cubic' | 'nearest'
%   Nmin       : minimum output size (default 100)
% Outputs:
%   Uq, xq, yq : interpolated field and axes ([1 x Nout] each)
%   Info       : InputSize, OutputSize, Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4 || isempty(method), method = 'spline'; end
if nargin < 5 || isempty(Nmin), Nmin = 100; end
[nB, nA, nT] = size(U);
Nout = max([nA, nB, Nmin]);
xq = linspace(x_mm(1), x_mm(end), Nout);
yq = linspace(y_mm(1), y_mm(end), Nout);
[Xq, Yq] = meshgrid(xq, yq);

valid = all(isfinite(U), 3);
Fm = griddedInterpolant({y_mm(:), x_mm(:)}, double(valid), 'linear', 'none');
vq = Fm(Yq, Xq) > 0.999;

Uq = nan(Nout, Nout, nT, 'like', U);
for t = 1:nT
    u = U(:, :, t);
    u(~valid) = mean(u(valid));                 % neutral fill (masked afterwards)
    if isreal(u)
        F = griddedInterpolant({y_mm(:), x_mm(:)}, u, method, 'none');
        q = F(Yq, Xq);
    else
        Fr = griddedInterpolant({y_mm(:), x_mm(:)}, real(u), method, 'none');
        Fi = griddedInterpolant({y_mm(:), x_mm(:)}, imag(u), method, 'none');
        q = complex(Fr(Yq, Xq), Fi(Yq, Xq));
    end
    q(~vq) = NaN;
    Uq(:, :, t) = q;
end
Info = struct('InputSize', [nB nA nT], 'OutputSize', [Nout Nout nT], 'Method', method);
end
