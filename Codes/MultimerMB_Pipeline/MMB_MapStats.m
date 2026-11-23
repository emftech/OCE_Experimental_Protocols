function Stats = MMB_MapStats(Maps, c_ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_MapStats
%
% Global statistics of the four R^2-masked elastograms of one acquisition.
% The phantom is homogeneous, so the ROI is the whole quality-masked map.
%
%   Ep = sigma_ROI / mu_ROI x 100 %          (thesis Eq. 29)
%   Ea = |mu_ROI - c_ref| / c_ref x 100 %    (thesis Eq. 28)
%   Ea_pix = mean_pixels |c - c_ref| / c_ref x 100 % (local, no cancellation)
%
% Inputs:
%   Maps  : struct with fields PG1D, PG1D_TV, PG2D, PG2D_TV
%   c_ref : reference speed [m/s] (optional; NaN -> Ea not computed)
% Output:
%   Stats : table, one row per estimator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    c_ref = NaN;
end

names = {'PG1D', 'PG1D_TV', 'PG2D', 'PG2D_TV'};
Method = {'1D PG'; '1D PG-TV'; '2D PG'; '2D PG-TV'};
n = numel(names);
[Median, Q1, Q3, Mean, Std, Ep, Ea, Ea_pix, ValidPx] = deal(nan(n,1));

for i = 1:n
    v = Maps.(names{i});
    v = v(isfinite(v));
    ValidPx(i) = numel(v);
    if isempty(v), continue; end
    Median(i) = median(v);
    Q1(i) = prctile(v, 25);
    Q3(i) = prctile(v, 75);
    Mean(i) = mean(v);
    Std(i) = std(v);
    Ep(i) = 100*Std(i)/Mean(i);
    Ea(i) = 100*abs(Mean(i) - c_ref)/c_ref;
    Ea_pix(i) = 100*mean(abs(v - c_ref))/c_ref;
end

Stats = table(Method, Median, Q1, Q3, Mean, Std, Ep, Ea, Ea_pix, ValidPx);

end
