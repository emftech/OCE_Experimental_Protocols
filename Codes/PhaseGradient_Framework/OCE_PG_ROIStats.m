function Stats = OCE_PG_ROIStats(Results, Inc_Mask, Bkg_Mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_PG_ROIStats
%
% ROI statistics of the masked elastograms in Results.Maps, using the same
% metric definitions as the simulation study (thesis Chapters 4-6):
%
%   Precision error : Ep  = sigma_ROI / mu_ROI * 100 %
%   Contrast-to-noise: CNR = |mu_bkg - mu_inc| / sqrt((sigma_bkg^2 + sigma_inc^2)/2)
%
% Accuracy error is not computed because experimental data has no ground
% truth; compare the ROI means against an independent reference instead
% (e.g. mechanical testing or group-speed estimates).
%
% Inputs:
%   Results  : output of OCE_PG_ProcessMB
%   Inc_Mask : logical [Ny x Nx] inclusion ROI (can be empty)
%   Bkg_Mask : logical [Ny x Nx] background ROI (can be empty)
%
% Output:
%   Stats : table, one row per method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = fieldnames(Results.Maps);
n = numel(names);

Mean_Inc = nan(n,1); Std_Inc = nan(n,1); Ep_Inc = nan(n,1); N_Inc = zeros(n,1);
Mean_Bkg = nan(n,1); Std_Bkg = nan(n,1); Ep_Bkg = nan(n,1); N_Bkg = zeros(n,1);
CNR = nan(n,1);

for i = 1:n

    E = Results.Maps.(names{i});

    if ~isempty(Inc_Mask)
        v = E(logical(Inc_Mask) & isfinite(E));
        N_Inc(i) = numel(v);
        Mean_Inc(i) = mean(v);
        Std_Inc(i) = std(v);
        Ep_Inc(i) = 100*Std_Inc(i)/Mean_Inc(i);
    end

    if ~isempty(Bkg_Mask)
        v = E(logical(Bkg_Mask) & isfinite(E));
        N_Bkg(i) = numel(v);
        Mean_Bkg(i) = mean(v);
        Std_Bkg(i) = std(v);
        Ep_Bkg(i) = 100*Std_Bkg(i)/Mean_Bkg(i);
    end

    CNR(i) = abs(Mean_Bkg(i) - Mean_Inc(i)) / sqrt((Std_Bkg(i)^2 + Std_Inc(i)^2)/2);
end

Stats = table(Mean_Inc, Std_Inc, Ep_Inc, N_Inc, Mean_Bkg, Std_Bkg, Ep_Bkg, N_Bkg, CNR, ...
    'RowNames', names);

end
