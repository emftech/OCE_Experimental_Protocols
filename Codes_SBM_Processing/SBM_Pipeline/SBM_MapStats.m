function Stats = SBM_MapStats(Maps, c_ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_MapStats
%
% Global statistics of the R^2-masked SBM elastograms (2D PG, 2D PG-TV),
% same definitions as MMB_MapStats (homogeneous phantom: ROI = whole map):
%   Ep = sigma/mu x 100 %               (thesis Eq. 29)
%   Ea = |mu - c_ref|/c_ref x 100 %     (thesis Eq. 28)
%   Ea_pix = mean |c - c_ref|/c_ref x 100 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2, c_ref = NaN; end
names = {'PG2D', 'PG2D_TV'};
Method = {'2D PG'; '2D PG-TV'};
n = numel(names);
[Median, Q1, Q3, Mean, Std, Ep, Ea, Ea_pix, ValidPx] = deal(nan(n, 1));
for i = 1:n
    v = Maps.(names{i});
    v = v(isfinite(v));
    ValidPx(i) = numel(v);
    if isempty(v), continue; end
    Median(i) = median(v);  Q1(i) = prctile(v, 25);  Q3(i) = prctile(v, 75);
    Mean(i) = mean(v);  Std(i) = std(v);
    Ep(i) = 100*Std(i)/Mean(i);
    Ea(i) = 100*abs(Mean(i) - c_ref)/c_ref;
    Ea_pix(i) = 100*mean(abs(v - c_ref))/c_ref;
end
Stats = table(Method, Median, Q1, Q3, Mean, Std, Ep, Ea, Ea_pix, ValidPx);
end
