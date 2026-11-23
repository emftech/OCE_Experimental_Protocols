%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% Main_Processing_Batch_MultimerMB_EMF
%
% Batch processing of the Multimeridian MB-mode experiments listed in
% Experiment_Log_06_10_2026.xlsx (phantom PC, QH 3-cycle, 1 Vpp, 12 x 12 mm):
%
%   * 2-meridian tests          : 220 / 285 / 500 us
%   * Multiple meridians 285 us : 4, 8, 16, 32 meridians x r1-r3 (thesis replica)
%   * Multiple meridians 500 us : 4, 8, 16, 32 meridians x r1-r3
%   * Frequency sweep           : 220-340 us (step 20 us), 8 meridians, r1
%
% Stage 1 - every experiment is processed with
%   Main_Processing_SingleExp_MultimerMB_EMF.m (raw data -> elastograms,
%   without and with directional filtering). Incomplete .bin files are
%   detected and skipped. The list can be split between several MATLAB
%   processes (WorkerID / NumWorkers) to use the 48 GB of RAM in parallel.
%
% Stage 2 - summary (Results/Summary_MMB_07-10-2026/):
%   * summary table of all experiments, estimators and filter options
%   * best repetition per (period, meridians) = the one with the largest
%     number of well-fitted pixels (R^2 >= 0.7) over the four estimators
%   * 4 x 4 elastogram grids (rows 4/8/16/32 meridians, columns 1D PG,
%     1D PG-TV, 2D PG, 2D PG-TV) for 285 us and 500 us, without/with DF
%   * accuracy (Ea) and precision (Ep) errors, thesis definitions
%   * dispersion analysis c(f) with 9 frequencies (7 sweep + 285/500 us
%     8-meridian runs), slope dc/df and Kelvin-Voigt fit
%   * spectra / time profiles of the sweep, SBM frames vs frequency map
%
% Workspace options (set before running):
%   OnlyExperiments = {...};   Overwrite = true;   SummaryOnly = true;
%   WorkerID = 1; NumWorkers = 2;   SkipSummary = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;

%% Configuration

B = struct();
B.RootDir    = fileparts(mfilename('fullpath'));   % repository root
B.LogFile    = fullfile(B.RootDir, 'Experiment_Log_06_10_2026.xlsx');
B.LogSheet   = 'PH_1150_503_06_10_2026';
B.DataDir    = fullfile(B.RootDir, 'Data');
B.ResultsDir = fullfile(B.RootDir, 'Results');
B.SummaryDir = fullfile(B.ResultsDir, 'Summary_MMB_07-10-2026');
B.CodesDir   = fullfile(B.RootDir, 'Codes');
B.SimCodesDir = fullfile(B.RootDir, '..', 'SBM_Protocol_OCE', 'Codes');   % simulator repo cloned next to this one

% Reference phase speed for the accuracy error (thesis Eq. 28).
% Primary: Rayleigh speed of the PC phantom from the quasi-static
% compression tests (Compute_Accuracy_Theory.m, 05-10-2026):
%   E0 = 0.365 MPa, rho = 1000 kg/m^3, nu -> 0.5: c_R = 0.9553 sqrt(E0/(3 rho))
B.Ref.Theory_cR = 0.9553*sqrt(0.365e6/(3*1000));      % 10.54 m/s
% Secondary: mean speed measured on the same phantom on 26-09-2026
% (all estimators, >= 4 meridians), rescaled from the 12 mm assumed then
% to the 10 mm written in those .bin headers (see Processing 07-10-2026).
B.Ref.PrevSummary = fullfile(B.RootDir, 'Previous', '26-09-2026', 'Results', ...
    'Summary_PG_Batch', 'PG_Summary.csv');
B.Ref.PrevScale = 10/12;

B.SweepPeriods_us = [285 500 220 240 260 280 300 320 340];   % SBM figure / dispersion

if ~exist('OnlyExperiments', 'var'), OnlyExperiments = {}; end
if ~exist('Overwrite', 'var'),       Overwrite = false;     end
if ~exist('SummaryOnly', 'var'),     SummaryOnly = false;   end
if ~exist('SkipSummary', 'var'),     SkipSummary = false;   end
if ~exist('WorkerID', 'var'),        WorkerID = 1;          end
if ~exist('NumWorkers', 'var'),      NumWorkers = 1;        end

addpath(genpath(B.CodesDir));
set(0, 'DefaultFigureVisible', 'off');

%% Experiment list from the log

Exps = readLog(B);
fprintf('%d experiments in the log.\n', numel(Exps));

%% Stage 1 - processing

if ~SummaryOnly
    todo = 1:numel(Exps);
    if ~isempty(OnlyExperiments)
        todo = find(ismember({Exps.Name}, OnlyExperiments));
    end
    % Largest files first, dealt round-robin to the workers
    [~, o] = sort([Exps(todo).Meridians], 'descend');
    todo = todo(o);
    todo = todo(WorkerID:NumWorkers:end);

    for e = todo
        E = Exps(e);
        matFile = fullfile(B.ResultsDir, E.Name, 'Results_MMB.mat');
        if ~E.FileOK
            fprintf('[%s] %s: %s -> skipped\n', datestr(now, 'HH:MM'), E.Name, E.FileNote); %#ok<TNOW1,DATST>
            continue;
        end
        if exist(matFile, 'file') && ~Overwrite
            fprintf('[%s] %s: already processed\n', datestr(now, 'HH:MM'), E.Name); %#ok<TNOW1,DATST>
            continue;
        end
        fprintf('[%s] ', datestr(now, 'HH:MM')); %#ok<TNOW1,DATST>
        try
            runSingle(E.Name);
        catch ME
            fprintf(2, '    FAILED: %s\n', ME.message);
            for k = 1:numel(ME.stack)
                fprintf(2, '      at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
            end
        end
        close all;
    end
end

%% Stage 2 - summary

if ~SkipSummary
    buildSummary(Exps, B);
end

%% ======================================================================= %
%                            LOCAL FUNCTIONS
% ======================================================================== %

function runSingle(ExpName)
% Runs the single-experiment script in this function's own workspace.
BatchMode = true; %#ok<NASGU>
Main_Processing_SingleExp_MultimerMB_EMF;
end

%% ----------------------------------------------------------------------- %
function Exps = readLog(B)
% Experiments of the 06-10-2026 log: name, cluster (group) and file check.

C = readcell(B.LogFile, 'Sheet', B.LogSheet);
nameCol = size(C, 2);
Exps = struct([]);
group = '';

for r = 1:size(C, 1)
    c1 = C{r, 1};
    if (ischar(c1) || isstring(c1)) && strlength(strtrim(string(c1))) > 0
        group = char(strtrim(string(c1)));
    end
    nm = C{r, nameCol};
    if ~(ischar(nm) || isstring(nm)) || ~contains(string(nm), '_OCE_')
        continue;
    end
    E = MMB_ParseExperimentName(char(nm));
    if contains(group, 'Tests')
        E.Group = '2-meridian tests';
    elseif contains(group, '285')
        E.Group = 'Multimeridian 285 us';
    elseif contains(group, '500')
        E.Group = 'Multimeridian 500 us';
    elseif contains(group, 'Frequency')
        E.Group = 'Frequency sweep';
    else
        E.Group = group;
    end
    E.LogGroup = group;

    hit = dir(fullfile(B.DataDir, '**', [E.Name '.bin']));
    if isempty(hit)
        E.FileOK = false;
        E.FileNote = 'file not found';
    else
        [~, ~, chk] = MMB_ReadBinHeader(fullfile(hit(1).folder, hit(1).name));
        E.FileOK = chk.IsComplete;
        if chk.IsComplete
            E.FileNote = 'complete';
        else
            E.FileNote = sprintf('incomplete file (%.2f of %.2f GB, %d of %d meridians)', ...
                chk.ActualBytes/1e9, chk.ExpectedBytes/1e9, chk.CompleteMeridians, E.Meridians);
        end
    end
    Exps = [Exps, E]; %#ok<AGROW>
end

end

%% ----------------------------------------------------------------------- %
function buildSummary(Exps, B)

S = B.SummaryDir;
if ~exist(S, 'dir'), mkdir(S); end
names = {'PG1D', 'PG1D_TV', 'PG2D', 'PG2D_TV'};
methods = {'1D PG', '1D PG-TV', '2D PG', '2D PG-TV'};
dfTags = {'noDF', 'DF'};
dfLabels = {'No directional filter', 'Directional filter (outward)'};

%% ---- References for the accuracy error
refPrev = previousReference(B);         % table: f_exc_Hz, cPrev
cTheory = B.Ref.Theory_cR;

%% ---- Load all results
R = struct();
T = table();
for e = 1:numel(Exps)
    E = Exps(e);
    f = fullfile(B.ResultsDir, E.Name, 'Results_MMB.mat');
    if ~exist(f, 'file'), continue; end
    L = load(f, 'Res', 'ResDF', 'Info', 'DFInfo');
    key = matlab.lang.makeValidName(E.Name);
    R.(key).Maps = {L.Res.Maps, L.ResDF.Maps};
    R.(key).x = L.Res.x_mm; R.(key).y = L.Res.y_mm;
    R.(key).Info = L.Info;
    R.(key).E = E;

    iPrev = find(abs(refPrev.f_exc_Hz - E.f_exc) < 1, 1);
    cPrev = NaN;
    if ~isempty(iPrev), cPrev = refPrev.cPrev(iPrev); end

    for d = 1:2
        St  = MMB_MapStats(R.(key).Maps{d}, cTheory);
        St2 = MMB_MapStats(R.(key).Maps{d}, cPrev);
        for i = 1:4
            row = table(string(E.Name), string(E.Group), E.Period_us, E.f_exc, L.Info.f_an, ...
                L.Info.Band.f_peak, L.Info.Band.BW_Hz, E.Meridians, string(E.Repetition), ...
                string(dfTags{d}), string(methods{i}), St.Median(i), St.Mean(i), St.Std(i), ...
                St.Ep(i), St.ValidPx(i), cTheory, St.Ea(i), St.Ea_pix(i), cPrev, St2.Ea(i), ...
                L.Info.Source_mm(1), L.Info.Source_mm(2), L.Info.Source.Fit_mm(1), L.Info.Source.Fit_mm(2), ...
                L.Info.Source.Coherence, L.Info.Source.Coherence_Centered, L.Info.Source.c_mps, ...
                median(L.DFInfo.OutwardEnergyFraction), ...
                'VariableNames', {'Experiment', 'Group', 'Period_us', 'f_exc_Hz', 'f_an_Hz', ...
                'f_peak_Hz', 'BW_Hz', 'Meridians', 'Repetition', 'Filter', 'Method', 'Median_mps', ...
                'Mean_mps', 'Std_mps', 'Ep_pct', 'ValidPx', 'cRef_Theory_mps', 'Ea_Theory_pct', ...
                'EaPix_Theory_pct', 'cRef_Prev_mps', 'Ea_Prev_pct', 'Source_x_mm', 'Source_y_mm', ...
                'SourceFit_x_mm', 'SourceFit_y_mm', ...
                'Source_Coherence', 'Source_CoherenceCentred', 'c_CylFit_mps', 'OutwardEnergyFraction'});
            T = [T; row]; %#ok<AGROW>
        end
    end
end
writetable(T, fullfile(S, 'MMB_Summary_All.csv'));
writetable(T, fullfile(S, 'MMB_Summary_All.xlsx'));

ExpT = struct2table(rmfield(Exps, {'RepIdx'}));
writetable(ExpT, fullfile(S, 'Experiment_List.csv'));

%% ---- Best repetition per (period, meridians)
Best = table();
for per = [285 500]
    for mer = [4 8 16 32]
        sel = T(T.Period_us == per & T.Meridians == mer & T.Filter == "noDF", :);
        if isempty(sel), continue; end
        reps = unique(sel.Repetition);
        score = arrayfun(@(r) sum(sel.ValidPx(sel.Repetition == r)), reps);
        [~, ib] = max(score);
        nm = sel.Experiment(find(sel.Repetition == reps(ib), 1));
        Best = [Best; table(per, mer, reps(ib), nm, score(ib), strjoin(reps + ":" + string(score), ', '), ...
            'VariableNames', {'Period_us', 'Meridians', 'BestRep', 'Experiment', 'R2Pixels_4methods', ...
            'AllReps'})]; %#ok<AGROW>
    end
end
writetable(Best, fullfile(S, 'Best_Repetition.csv'));
disp(Best);

%% ---- 4 x 4 elastogram grids (best repetition), without / with DF
for per = [285 500]
    rows = Best(Best.Period_us == per, :);
    if isempty(rows), continue; end
    cells = cell(height(rows), 4, 2);
    srcs = cell(height(rows), 1);
    labels = cell(height(rows), 1);
    for r = 1:height(rows)
        key = matlab.lang.makeValidName(rows.Experiment(r));
        for d = 1:2
            for i = 1:4
                cells{r, i, d} = R.(key).Maps{d}.(names{i});
            end
        end
        srcs{r} = R.(key).Info.Source_mm;
        labels{r} = sprintf('%d (%s)', rows.Meridians(r), rows.BestRep(r));
        x = R.(key).x; y = R.(key).y;
    end
    CL = MMB_SpeedCLim(cells(:));
    for d = 1:2
        fig = MMB_PlotElastogramGrid(cells(:,:,d), x, y, labels, methods, CL, ...
            sprintf('Multimeridian MB Mode: Elastogram Comparison - f_{exc} = %.0f Hz, %s', ...
            1e6/per, dfLabels{d}), 'Source_mm', srcs);
        MMB_SaveFigure(fig, fullfile(S, sprintf('Grid4x4_%dus_%s', per, dfTags{d})));
        close(fig);
    end
end

%% ---- Accuracy and precision vs number of meridians (mean +- SD over reps)
metricNames = {'Ea_Theory_pct', 'Ep_pct', 'Ea_Prev_pct', 'Median_mps'};
metricLabels = {'E_a [%] (ref. c_R theory)', 'E_p [%]', 'E_a [%] (ref. 26-09 mean)', 'Median phase speed [m/s]'};
mers = [4 8 16 32];
Agg = table();
for per = [285 500]
    for d = 1:2
        for i = 1:4
            for mer = mers
                sel = T(T.Period_us == per & T.Meridians == mer & T.Filter == dfTags{d} & T.Method == methods{i}, :);
                if isempty(sel), continue; end
                row = table(per, string(dfTags{d}), string(methods{i}), mer, height(sel), ...
                    'VariableNames', {'Period_us', 'Filter', 'Method', 'Meridians', 'nReps'});
                for q = 1:numel(metricNames)
                    v = sel.(metricNames{q});
                    row.([metricNames{q} '_mean']) = mean(v, 'omitnan');
                    row.([metricNames{q} '_sd']) = std(v, 'omitnan');
                end
                Agg = [Agg; row]; %#ok<AGROW>
            end
        end
    end
end
writetable(Agg, fullfile(S, 'Metrics_vs_Meridians.csv'));
writetable(Agg, fullfile(S, 'Metrics_vs_Meridians.xlsx'));

cols = lines(4);
for q = 1:numel(metricNames)
    for per = [285 500]
        fig = figure('Color', 'w', 'Position', [40 40 1300 520]);
        tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
        for d = 1:2
            ax = nexttile(tl); hold(ax, 'on');
            for i = 1:4
                sel = Agg(Agg.Period_us == per & Agg.Filter == dfTags{d} & Agg.Method == methods{i}, :);
                if isempty(sel), continue; end
                errorbar(ax, sel.Meridians, sel.([metricNames{q} '_mean']), sel.([metricNames{q} '_sd']), ...
                    '-o', 'Color', cols(i,:), 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'w', ...
                    'CapSize', 8, 'DisplayName', methods{i});
            end
            if strcmp(metricNames{q}, 'Median_mps')
                yline(ax, cTheory, 'k--', sprintf('c_R theory = %.2f m/s', cTheory), 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');
            end
            hold(ax, 'off'); grid(ax, 'on'); grid(ax, 'minor'); box(ax, 'on');
            set(ax, 'XTick', mers, 'FontSize', 13, 'FontWeight', 'bold');
            xlim(ax, [2 34]);
            xlabel(ax, 'Number of Meridians'); ylabel(ax, metricLabels{q});
            title(ax, dfLabels{d});
            legend(ax, 'Location', 'best');
        end
        linkaxes(tl.Children(end:-1:1), 'y');
        title(tl, sprintf('Multimeridian MB Mode - %s - f_{exc} = %.0f Hz (mean \\pm SD over repetitions)', ...
            strrep(metricLabels{q}, ' [%]', ''), 1e6/per), 'FontSize', 15, 'FontWeight', 'bold');
        MMB_SaveFigure(fig, fullfile(S, sprintf('%s_vs_Meridians_%dus', metricNames{q}, per)));
        close(fig);
    end
end

%% ---- Dispersion analysis (9 frequencies, 8 meridians)
Disp = table();
for per = B.SweepPeriods_us
    for d = 1:2
        for i = 1:4
            sel = T(T.Period_us == per & T.Meridians == 8 & T.Filter == dfTags{d} & T.Method == methods{i}, :);
            if isempty(sel), continue; end
            Disp = [Disp; table(per, mean(sel.f_an_Hz), string(dfTags{d}), string(methods{i}), height(sel), ...
                mean(sel.Median_mps), std(sel.Median_mps), mean(sel.Ep_pct), ...
                'VariableNames', {'Period_us', 'f_Hz', 'Filter', 'Method', 'nReps', 'c_mps', 'c_sd_mps', ...
                'Ep_pct'})]; %#ok<AGROW>
        end
    end
end
Disp = sortrows(Disp, {'Filter', 'Method', 'f_Hz'});

% Linear fit c = c0 + (dc/df) f and Kelvin-Voigt Rayleigh-wave fit
rho = 1000;
kvModel = @(p, w) 0.9553*sqrt(2*(p(1)^2 + (w*p(2)).^2) ./ (rho*(p(1) + sqrt(p(1)^2 + (w*p(2)).^2))));
% c(f) is not monotonic in this phantom (rises up to ~3.5 kHz, then
% decreases), so the slope dc/df is fitted over the full range and over
% the low (2.0-3.6 kHz) and high (3.5-4.6 kHz) sub-bands.
Ranges = struct('Name', {'All 2.0-4.6 kHz', 'Low 2.0-3.6 kHz', 'High 3.5-4.6 kHz'}, 'Lim', {[0 5e3], [0 3.6e3], [3.5e3 5e3]});
Fit = table();
for d = 1:2
    for i = 1:4
        for g = 1:numel(Ranges)
            sel = Disp(Disp.Filter == dfTags{d} & Disp.Method == methods{i} & Disp.f_Hz >= Ranges(g).Lim(1) & Disp.f_Hz <= Ranges(g).Lim(2), :);
            if height(sel) < 3, continue; end
            f = sel.f_Hz; c = sel.c_mps;
            mdl = fitlm(f/1e3, c);
            ci = coefCI(mdl);
            p = [NaN NaN];
            if g == 1
                w = 2*pi*f;
                p = abs(fminsearch(@(p) sum((kvModel(abs(p), w) - c).^2), [median(c)^2*rho 0.5]));
            end
            Fit = [Fit; table(string(dfTags{d}), string(methods{i}), string(Ranges(g).Name), height(sel), ...
                mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2), ci(2,1), ci(2,2), ...
                mdl.Rsquared.Ordinary, mdl.Coefficients.pValue(2), p(1), p(2), 3*p(1)/1e6, ...
                'VariableNames', {'Filter', 'Method', 'Range', 'nPoints', 'c0_mps', 'dcdf_mps_per_kHz', ...
                'dcdf_CI95_lo', 'dcdf_CI95_hi', 'R2', 'p_slope', 'mu_KV_Pa', 'eta_KV_Pas', 'E_KV_MPa'})]; %#ok<AGROW>
        end
    end
end
writetable(Disp, fullfile(S, 'Dispersion_Points.csv'));
writetable(Fit, fullfile(S, 'Dispersion_Fits.csv'));
writetable(Disp, fullfile(S, 'Dispersion.xlsx'), 'Sheet', 'Points');
writetable(Fit, fullfile(S, 'Dispersion.xlsx'), 'Sheet', 'Fits');
disp(Fit);

fig = figure('Color', 'w', 'Position', [40 40 1400 560]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
ff = linspace(1.8e3, 4.8e3, 200);
for d = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    for i = 1:4
        sel = Disp(Disp.Filter == dfTags{d} & Disp.Method == methods{i}, :);
        if isempty(sel), continue; end
        isMM = ismember(sel.Period_us, [285 500]);
        errorbar(ax, sel.f_Hz/1e3, sel.c_mps, sel.c_sd_mps, 'o', 'Color', cols(i,:), ...
            'MarkerFaceColor', cols(i,:), 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', methods{i});
        plot(ax, sel.f_Hz(isMM)/1e3, sel.c_mps(isMM), 's', 'Color', 'k', 'MarkerSize', 13, ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');
        for g = 2:3
            fr = Fit(Fit.Filter == dfTags{d} & Fit.Method == methods{i} & Fit.Range == Ranges(g).Name, :);
            if isempty(fr), continue; end
            fg = linspace(max(Ranges(g).Lim(1), 1.95e3), min(Ranges(g).Lim(2), 4.6e3), 50);
            ls = {'-', '--'};
            plot(ax, fg/1e3, fr.c0_mps + fr.dcdf_mps_per_kHz*fg/1e3, ls{g-1}, 'Color', cols(i,:), 'LineWidth', 1.4, ...
                'HandleVisibility', 'off');
        end
    end
    yline(ax, cTheory, 'k--', 'c_R (quasi-static compression)', 'HandleVisibility', 'off', 'LineWidth', 1.2);
    hold(ax, 'off'); grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Frequency [kHz]'); ylabel(ax, 'Phase speed [m/s]');
    title(ax, dfLabels{d});
    set(ax, 'FontSize', 13, 'FontWeight', 'bold');
    legend(ax, 'Location', 'northwest', 'FontSize', 9);
end
linkaxes(tl.Children(end:-1:1), 'y');
title(tl, {'Dispersion - 8 meridians, 9 excitation frequencies', ...
    '(squares: 285/500 us runs, mean \pm SD over 3 reps; solid: fit 2.0-3.6 kHz, dashed: fit 3.5-4.6 kHz)'}, 'FontSize', 14, 'FontWeight', 'bold');
MMB_SaveFigure(fig, fullfile(S, 'Dispersion_c_vs_f'));
close(fig);

%% ---- Frequency sweep: spectra, time profiles, 2D PG-TV elastograms
sweepNames = strings(0); sweepF = [];
for per = sort(B.SweepPeriods_us)
    sel = T(T.Period_us == per & T.Meridians == 8 & T.Repetition == "r1" & T.Filter == "noDF" & T.Method == "2D PG-TV", :);
    if isempty(sel), continue; end
    sweepNames(end+1) = sel.Experiment(1); %#ok<AGROW>
    sweepF(end+1) = sel.f_exc_Hz(1); %#ok<AGROW>
end
[sweepF, o] = sort(sweepF); sweepNames = sweepNames(o);
nSw = numel(sweepNames);
nc = 3; nr = ceil(nSw/nc);

fig = figure('Color', 'w', 'Position', [40 40 1500 950]);
tl = tiledlayout(nr, nc, 'TileSpacing', 'compact');
for k = 1:nSw
    key = matlab.lang.makeValidName(sweepNames(k)); I = R.(key).Info;
    ax = nexttile(tl);
    plot(ax, I.Spectrum_f/1e3, I.Spectrum_PdB, 'k', 'LineWidth', 1.3); hold(ax, 'on');
    patch(ax, [I.Band.f_lo I.Band.f_hi I.Band.f_hi I.Band.f_lo]/1e3, [-40 -40 3 3], [0.3 0.75 0.93], ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');
    xline(ax, I.f_exc/1e3, 'r--', 'LineWidth', 1.2);
    hold(ax, 'off'); grid(ax, 'on'); xlim(ax, [0 8]); ylim(ax, [-40 3]);
    title(ax, sprintf('%.0f Hz (%g \\mus): peak %.0f Hz, BW_{-3dB} %.0f Hz', I.f_exc, 1e6/I.f_exc, ...
        I.Band.f_peak, I.Band.BW3dB_Hz), 'FontSize', 10);
    if k > nSw - nc, xlabel(ax, 'Frequency [kHz]'); end
    if mod(k-1, nc) == 0, ylabel(ax, 'Power [dB]'); end
end
title(tl, 'Mean surface spectra (8 meridians, r1) and narrowband pass bands', 'FontWeight', 'bold');
MMB_SaveFigure(fig, fullfile(S, 'Sweep_Spectra'));
close(fig);

fig = figure('Color', 'w', 'Position', [40 40 1500 950]);
tl = tiledlayout(nr, nc, 'TileSpacing', 'compact');
for k = 1:nSw
    L = load(fullfile(B.ResultsDir, sweepNames(k), 'Results_MMB.mat'), 'loaded_phases_Border', 'Frames', 'Time', 'Info');
    p = L.Info.RepresentativePos;
    ax = nexttile(tl);
    plot(ax, L.Time, L.loaded_phases_Border(p,:), 'Color', [0.4 0.4 0.4]); hold(ax, 'on');
    plot(ax, L.Time, L.Frames(p,:), 'r', 'LineWidth', 1.4); hold(ax, 'off'); grid(ax, 'on');
    xlim(ax, [L.Time(1) L.Time(end)]);
    title(ax, sprintf('%.0f Hz', L.Info.f_exc), 'FontSize', 11);
    if k > nSw - nc, xlabel(ax, 'time (ms)'); end
    if mod(k-1, nc) == 0, ylabel(ax, 'v_z [mm/s]'); end
end
title(tl, 'Surface particle velocity: raw (grey) and narrowband (red)', 'FontWeight', 'bold');
MMB_SaveFigure(fig, fullfile(S, 'Sweep_TimeProfiles'));
close(fig);

for d = 1:2
    for i = [2 4]
        cells = cell(nr, nc); srcs = cell(nr, 1);
        allm = {};
        for k = 1:nSw
            key = matlab.lang.makeValidName(sweepNames(k));
            cells{ceil(k/nc), mod(k-1, nc)+1} = R.(key).Maps{d}.(names{i});
            allm{end+1} = R.(key).Maps{d}.(names{i}); %#ok<AGROW>
        end
        key1 = matlab.lang.makeValidName(sweepNames(1));
        fig = MMB_PlotElastogramGrid(cells, R.(key1).x, R.(key1).y, ...
            arrayfun(@(r) sprintf('%.0f-%.0f Hz', sweepF((r-1)*nc+1), sweepF(min(r*nc, nSw))), 1:nr, 'UniformOutput', false), ...
            arrayfun(@(c) sprintf('col %d', c), 1:nc, 'UniformOutput', false), MMB_SpeedCLim(allm), ...
            sprintf('Frequency sweep (8 meridians, r1) - %s, %s', methods{i}, dfLabels{d}), 'RowSuffix', '');
        % Put the frequency of each panel as its title
        axs = flipud(findobj(fig, 'Type', 'axes'));
        for k = 1:min(nSw, numel(axs))
            title(axs(k), sprintf('%.0f Hz', sweepF(k)), 'FontSize', 13, 'FontWeight', 'bold');
        end
        MMB_SaveFigure(fig, fullfile(S, sprintf('Sweep_Grid_%s_%s', names{i}, dfTags{d})));
        close(fig);
    end
end

%% ---- SBM: feasible synchronized frames vs frequency (9 periods)
fig = MMB_PlotSBMFrames(B.SweepPeriods_us, 5, fullfile(S, 'SBM_Frames_vs_Frequency_9periods'));
close(fig);
SBM = table();
for per = B.SweepPeriods_us
    M = per/5;
    SBM = [SBM; table(per, 1e6/per, M, string(strjoin(string(divisors(M)), ', ')), ...
        'VariableNames', {'Period_us', 'f_exc_Hz', 'M', 'FeasibleFrames'})]; %#ok<AGROW>
end
writetable(SBM, fullfile(S, 'SBM_Feasible_Frames.csv'));

save(fullfile(S, 'Summary_MMB.mat'), 'T', 'Best', 'Agg', 'Disp', 'Fit', 'SBM', 'refPrev', 'cTheory', 'B', 'Exps');
fprintf('Summary written to %s\n', S);

end

%% ----------------------------------------------------------------------- %
function ref = previousReference(B)
% Mean phase speed of the PC phantom measured on 26-09-2026 (four
% estimators, >= 4 meridians), rescaled to the 10 mm header scan length.
ref = table([2000; 1e6/285], [NaN; NaN], 'VariableNames', {'f_exc_Hz', 'cPrev'});
if ~exist(B.Ref.PrevSummary, 'file'), return; end
P = readtable(B.Ref.PrevSummary);
for k = 1:height(ref)
    sel = abs(P.f_exc_Hz - ref.f_exc_Hz(k)) < 1 & P.Meridians >= 4;
    ref.cPrev(k) = mean(P.Mean_mps(sel), 'omitnan')*B.Ref.PrevScale;
end
end
