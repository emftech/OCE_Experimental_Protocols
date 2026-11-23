%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% Main_Processing_Batch_SBM_EMF
%
% Batch processing of the Synchronized BM-mode (SBM) experiments of the Twin
% OCT (Experiment_Log_SBM_*.xlsx: phantom PC, Rigol 1-cycle bursts locked to
% the A-line clock, 1 Vpp, 12 x 12 mm), mirroring
% Main_Processing_Batch_MultimerMB_EMF.m:
%
%   * Tests at 285 us    : 37 x 21 B-scans, 57 frames, r1-r5 (repeatability)
%   * B-scan sweep 285 us: 11 / 21 / 41 / 81 B-scans (~ 4 / 8 / 16 / 32 meridians)
%   * Frequency sweep    : 220-340 us (20 us steps), 285 us and 500 us, 21 B-scans
%
% The experiment list is built from the .bin files present in
% SBM_Experiments/Data (names follow the log convention), so the batch does
% not depend on the Excel file (whose formula cells may have no cached
% values).
%
% Stage 1 - every experiment with Main_Processing_SingleExp_SBM_EMF.m
%   (raw data -> 2D PG / 2D PG-TV, without and with directional filtering).
%   Incomplete files are skipped. WorkerID / NumWorkers split the list
%   between several MATLAB processes, as in the MB batch.
% Stage 2 - summary (SBM_Experiments/Results/Summary_SBM/):
%   * SBM_Summary_All.csv/.xlsx (every experiment x estimator x filter)
%   * repeatability of the 285 us tests (mean +- SD, grids of all reps)
%   * B-scan sweep: elastogram grids 11/21/41/81 x estimators, metrics vs
%     N_B, and comparison with the Multimeridian MB metrics vs meridians
%     (4/8/16/32, Summary_MMB.mat) including acquisition time
%   * dispersion c(f) of the frequency sweep (+ MB points), dc/df fits
%   * joint (session) needle position from the summed coherence maps
%
% Workspace options (set before running):
%   OnlyExperiments = {...};  Overwrite = true;  SummaryOnly = true;
%   SkipSummary = true;  WorkerID = 1; NumWorkers = 2;
%   CfgOverride = struct(...)  (passed to every single-experiment run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('KeepWorkspace', 'var') || ~KeepWorkspace
    clc; close all;
end

%% Configuration

B = struct();
B.RootDir     = fileparts(fileparts(mfilename('fullpath')));   % repository root
B.ExpDir      = fullfile(B.RootDir, 'SBM_Experiments');
B.DataDir     = fullfile(B.ExpDir, 'Data');
B.ResultsDir  = fullfile(B.ExpDir, 'Results');
B.SummaryDir  = fullfile(B.ResultsDir, 'Summary_SBM');
B.CodesDir    = fullfile(B.RootDir, 'Codes');
B.SBMCodesDir = fullfile(B.RootDir, 'Codes_SBM_Processing');
B.MBSummary   = fullfile(B.RootDir, 'Results', 'Summary_MMB_07-10-2026', 'Summary_MMB.mat');
B.Ref.Theory_cR = 0.9553*sqrt(0.365e6/(3*1000));     % 10.54 m/s (same reference as MB)
B.TestPeriod_us = 285;  B.TestNB = 21;
B.MBAnalog = [11 4; 21 8; 41 16; 81 32];              % SBM B-scans <-> MB meridians (thesis Table 6)
B.MBAcqPerMeridian_s = 100*(450+50)*5e-6;             % 06-10 MB protocol: 100 positions x 500 A-lines

if ~exist('OnlyExperiments', 'var'), OnlyExperiments = {}; end
if ~exist('Overwrite', 'var'),       Overwrite = false;     end
if ~exist('SummaryOnly', 'var'),     SummaryOnly = false;   end
if ~exist('SkipSummary', 'var'),     SkipSummary = false;   end
if ~exist('WorkerID', 'var'),        WorkerID = 1;          end
if ~exist('NumWorkers', 'var'),      NumWorkers = 1;        end
if ~exist('CfgOverride', 'var'),     CfgOverride = struct(); end
if isfield(CfgOverride, 'DataDir'),    B.DataDir = CfgOverride.DataDir; end
if isfield(CfgOverride, 'ResultsDir'), B.ResultsDir = CfgOverride.ResultsDir; B.SummaryDir = fullfile(B.ResultsDir, 'Summary_SBM'); end

addpath(genpath(B.CodesDir));  addpath(genpath(B.SBMCodesDir));
set(0, 'DefaultFigureVisible', 'off');

%% Experiment list from the data folder

Exps = listExperiments(B);
fprintf('%d SBM acquisitions in %s.\n', numel(Exps), B.DataDir);

%% Stage 1 - processing

if ~SummaryOnly
    todo = 1:numel(Exps);
    if ~isempty(OnlyExperiments), todo = find(ismember({Exps.Name}, OnlyExperiments)); end
    [~, o] = sort([Exps(todo).Bytes], 'descend');      % largest files first
    todo = todo(o);
    todo = todo(WorkerID:NumWorkers:end);
    for e = todo
        E = Exps(e);
        matFile = fullfile(B.ResultsDir, E.Name, 'Results_SBM.mat');
        if ~E.FileOK
            fprintf('[%s] %s: %s -> skipped\n', char(datetime('now', 'Format', 'HH:mm')), E.Name, E.FileNote);
            continue;
        end
        if exist(matFile, 'file') && ~Overwrite
            fprintf('[%s] %s: already processed\n', char(datetime('now', 'Format', 'HH:mm')), E.Name);
            continue;
        end
        fprintf('[%s] ', char(datetime('now', 'Format', 'HH:mm')));
        try
            runSingle(E.Name, CfgOverride);
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

function runSingle(ExpName, CfgOverride) %#ok<INUSD>
BatchMode = true; %#ok<NASGU>
Main_Processing_SingleExp_SBM_EMF;
end

%% ----------------------------------------------------------------------- %
function Exps = listExperiments(B)
L = dir(fullfile(B.DataDir, '**', '*.bin'));
Exps = struct('Name', {}, 'File', {}, 'Bytes', {}, 'FileOK', {}, 'FileNote', {}, ...
    'Period_us', {}, 'NB', {}, 'NA', {}, 'NF', {}, 'Rep', {}, 'Group', {});
for i = 1:numel(L)
    [~, nm] = fileparts(L(i).name);
    try
        E = SBM_ParseExperimentName(nm);
    catch
        continue;                                      % not an SBM file
    end
    f = fullfile(L(i).folder, L(i).name);
    [~, ~, Ck] = SBM_ReadBinHeader(f);
    note = '';
    if ~Ck.IsComplete
        note = sprintf('incomplete (%d frames on disk)', Ck.CompleteFrames);
    end
    if E.Period_us == B.TestPeriod_us && E.NB == B.TestNB
        grp = 'Tests';
    elseif E.Period_us == B.TestPeriod_us
        grp = 'BscanSweep';
    else
        grp = 'FreqSweep';
    end
    Exps(end+1) = struct('Name', nm, 'File', f, 'Bytes', L(i).bytes, 'FileOK', Ck.IsComplete, ...
        'FileNote', note, 'Period_us', E.Period_us, 'NB', E.NB, 'NA', E.NA, 'NF', E.NF, ...
        'Rep', E.RepIdx, 'Group', grp); %#ok<AGROW>
end
end

%% ----------------------------------------------------------------------- %
function buildSummary(Exps, B)
S = B.SummaryDir;
if ~exist(S, 'dir'), mkdir(S); end
fprintf('\nSummary -> %s\n', S);

R = struct([]);
for e = 1:numel(Exps)
    f = fullfile(B.ResultsDir, Exps(e).Name, 'Results_SBM.mat');
    if ~exist(f, 'file'), continue; end
    L = load(f, 'Stats', 'StatsDF', 'Info', 'Res', 'ResDF');
    L.Exp = Exps(e);
    R = [R, L]; %#ok<AGROW>
end
if isempty(R), fprintf('No processed experiments yet.\n'); return; end

% ---- table of everything (MB column names where they exist)
filtersAll = ["noDF", "DF"];
rows = {};
for i = 1:numel(R)
    I = R(i).Info;  X = R(i).Exp;
    if isfield(I, 'Detection')
        det = {I.Detection.CohExc, I.Detection.Threshold, I.Detection.Detected};
    else
        det = {NaN, NaN, false};
    end
    for d = 1:2
        flt = filtersAll(d);
        St = SBM_MapStats(ifelse(d == 1, R(i).Res.Maps, R(i).ResDF.Maps), B.Ref.Theory_cR);
        Rn = 1.2;  Mg = 1;                                  % far-field ROI r >= R + margin
        if isfield(I.Cfg, 'NeedleRadius_mm'), Rn = I.Cfg.NeedleRadius_mm; end
        if isfield(I.Cfg, 'FarFieldMargin_mm'), Mg = I.Cfg.FarFieldMargin_mm; end
        Sf = SBM_FarFieldStats(ifelse(d == 1, R(i).Res.Maps, R(i).ResDF.Maps), R(i).Res.x_mm, ...
            R(i).Res.y_mm, I.Source_mm, Rn, Mg, 0.5, B.Ref.Theory_cR);
        for k = 1:height(St)
            rows(end+1, :) = {string(X.Name), string(X.Group), X.Period_us, I.f_exc, X.NA, X.NB, X.NF, ...
                X.Rep, flt, string(St.Method{k}), St.Median(k), St.Mean(k), St.Std(k), St.Ep(k), ...
                St.ValidPx(k), B.Ref.Theory_cR, St.Ea(k), St.Ea_pix(k), I.Spectrum.SNR_dB, ...
                det{:}, I.Source_mm(1), I.Source_mm(2), I.Source.Coherence, I.AcqTime_s, ...
                Rn + Mg, Sf.Median(k), Sf.Ep(k), Sf.Ea(k), Sf.ValidPx(k)}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'Experiment', 'Group', 'Period_us', 'f_exc_Hz', 'NA', ...
    'NB', 'NF', 'Repetition', 'Filter', 'Method', 'Median_mps', 'Mean_mps', 'Std_mps', 'Ep_pct', ...
    'ValidPx', 'cRef_Theory_mps', 'Ea_Theory_pct', 'EaPix_Theory_pct', 'SNR_exc_dB', ...
    'ExcCoherence', 'ExcCohNullP95', 'ExcDetected', 'Source_x_mm', 'Source_y_mm', 'Source_Coherence', 'AcqTime_s', ...
    'FarROI_rmin_mm', 'Median_Far_mps', 'Ep_Far_pct', 'Ea_Far_Theory_pct', 'ValidPx_Far'});
writetable(T, fullfile(S, 'SBM_Summary_All.csv'));
writetable(T, fullfile(S, 'SBM_Summary_All.xlsx'));
nNo = numel(unique(T.Experiment(~T.ExcDetected)));
if nNo > 0
    warning('SBM:NoExcitation', ['%d of %d experiments show NO excitation-locked wave ', ...
        '(ExcDetected = false): their speeds come from noise (they are still plotted; do not interpret them).'], ...
        nNo, numel(unique(T.Experiment)));
end

methods = ["2D PG", "2D PG-TV"];  filters = ["noDF", "DF"];

% ---- repeatability of the 285 us tests
isTest = arrayfun(@(r) strcmp(r.Exp.Group, 'Tests'), R);
Rep = table();
if any(isTest)
    for fl = filters
        for me = methods
            sel = T.Group == "Tests" & T.Filter == fl & T.Method == me;
            Rep = [Rep; table(fl, me, nnz(sel), mean(T.Median_mps(sel)), std(T.Median_mps(sel)), ...
                100*std(T.Median_mps(sel))/mean(T.Median_mps(sel)), mean(T.Ea_Theory_pct(sel)), ...
                mean(T.Ep_pct(sel)), mean(T.ValidPx(sel)), 'VariableNames', {'Filter', 'Method', 'nReps', ...
                'Median_mps_mean', 'Median_mps_sd', 'CV_reps_pct', 'Ea_Theory_pct_mean', 'Ep_pct_mean', ...
                'ValidPx_mean'})]; %#ok<AGROW>
        end
    end
    writetable(Rep, fullfile(S, 'Tests_Repeatability.csv'));
    disp(Rep);
    it = find(isTest);  [~, o] = sort(arrayfun(@(r) r.Exp.Rep, R(it)));  it = it(o);
    for d = 1:2
        MapCell = cell(numel(it), 2);
        for i = 1:numel(it)
            Mp = ifelse(d == 1, R(it(i)).Res.Maps, R(it(i)).ResDF.Maps);
            MapCell(i, :) = {Mp.PG2D, Mp.PG2D_TV};
        end
        fig = MMB_PlotElastogramGrid(MapCell, R(it(1)).Res.x_mm, R(it(1)).Res.y_mm, ...
            arrayfun(@(r) sprintf('r%d', r.Exp.Rep), R(it), 'UniformOutput', false), ...
            {'2D PG', '2D PG-TV'}, MMB_SpeedCLim(MapCell), ...
            sprintf('SBM tests 285 us, 37 x 21 - repeatability (%s)', filters(d)), 'RowSuffix', '', ...
            'Source_mm', arrayfun(@(r) r.Info.Source_mm, R(it), 'UniformOutput', false));
        MMB_SaveFigure(fig, fullfile(S, sprintf('Tests_Grid_%s', filters(d))));
    end
end

% best test repetition = most well-fitted pixels over both estimators (no DF), as in MB
bestTest = [];
if any(isTest)
    it = find(isTest);
    px = arrayfun(@(r) sum(r.Stats.ValidPx), R(it));
    [~, k] = max(px);  bestTest = it(k);
end

% ---- B-scan sweep (285 us): 11 / 21 (best test) / 41 / 81
NBs = B.MBAnalog(:, 1)';
idxNB = nan(size(NBs));
for q = 1:numel(NBs)
    if NBs(q) == B.TestNB && ~isempty(bestTest)
        idxNB(q) = bestTest;
    else
        c = find(arrayfun(@(r) r.Exp.Period_us == B.TestPeriod_us && r.Exp.NB == NBs(q), R), 1);
        if ~isempty(c), idxNB(q) = c; end
    end
end
Agg = table();
if nnz(isfinite(idxNB)) >= 2
    have = find(isfinite(idxNB));
    for d = 1:2
        MapCell = cell(numel(have), 2);
        for q = 1:numel(have)
            Mp = ifelse(d == 1, R(idxNB(have(q))).Res.Maps, R(idxNB(have(q))).ResDF.Maps);
            MapCell(q, :) = {Mp.PG2D, Mp.PG2D_TV};
        end
        fig = MMB_PlotElastogramGrid(MapCell, R(idxNB(have(1))).Res.x_mm, R(idxNB(have(1))).Res.y_mm, ...
            arrayfun(@(n) sprintf('%d', n), NBs(have), 'UniformOutput', false), {'2D PG', '2D PG-TV'}, ...
            MMB_SpeedCLim(MapCell), sprintf('SBM B-scan sweep, 285 us (%s)', filters(d)), ...
            'RowSuffix', 'B-scans', 'Source_mm', arrayfun(@(q) R(idxNB(q)).Info.Source_mm, have, 'UniformOutput', false));
        MMB_SaveFigure(fig, fullfile(S, sprintf('BscanSweep_Grid_%s', filters(d))));
    end
    for fl = filters
        for me = methods
            for q = have
                X = R(idxNB(q)).Exp;
                if NBs(q) == B.TestNB      % mean +- SD over the test repetitions
                    sel = T.Group == "Tests" & T.Filter == fl & T.Method == me;
                else
                    sel = T.Experiment == string(X.Name) & T.Filter == fl & T.Method == me;
                end
                Agg = [Agg; table(B.TestPeriod_us, fl, me, NBs(q), B.MBAnalog(q, 2), nnz(sel), ...
                    mean(T.Ea_Theory_pct(sel)), std(T.Ea_Theory_pct(sel)), mean(T.Ep_pct(sel)), ...
                    std(T.Ep_pct(sel)), mean(T.Median_mps(sel)), std(T.Median_mps(sel)), ...
                    mean(T.ValidPx(sel)), mean(T.AcqTime_s(sel)), B.MBAnalog(q, 2)*B.MBAcqPerMeridian_s, ...
                    'VariableNames', {'Period_us', 'Filter', 'Method', 'NB', 'MB_Meridians', 'nReps', ...
                    'Ea_Theory_pct_mean', 'Ea_Theory_pct_sd', 'Ep_pct_mean', 'Ep_pct_sd', ...
                    'Median_mps_mean', 'Median_mps_sd', 'ValidPx_mean', 'AcqTime_s', 'MB_AcqTime_s'})]; %#ok<AGROW>
            end
        end
    end
    writetable(Agg, fullfile(S, 'Metrics_vs_Bscans.csv'));
    writetable(Agg, fullfile(S, 'Metrics_vs_Bscans.xlsx'));

    MB = [];
    if exist(B.MBSummary, 'file'), MB = load(B.MBSummary, 'Agg', 'Disp'); end
    metric = {'Ea_Theory_pct', 'Ep_pct', 'Median_mps'};
    mlab = {'E_a [%] (c_{ref} = 10.54 m/s)', 'E_p [%]', 'Median speed [m/s]'};
    fig = figure('Color', 'w', 'Position', [40 40 1500 900]);
    tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact');
    for d = 1:2
        for q = 1:3
            ax = nexttile(tl); hold(ax, 'on');
            for me = methods
                a = Agg(Agg.Filter == filters(d) & Agg.Method == me, :);
                errorbar(ax, a.MB_Meridians, a.([metric{q} '_mean']), a.([metric{q} '_sd']), '-o', ...
                    'LineWidth', 1.6, 'DisplayName', "SBM " + me);
                if ~isempty(MB)
                    b = MB.Agg(MB.Agg.Period_us == B.TestPeriod_us & MB.Agg.Filter == filters(d) & MB.Agg.Method == me, :);
                    errorbar(ax, b.Meridians, b.([metric{q} '_mean']), b.([metric{q} '_sd']), '--s', ...
                        'LineWidth', 1.2, 'DisplayName', "MB " + me);
                end
            end
            hold(ax, 'off'); grid(ax, 'on'); set(ax, 'XScale', 'log', 'XTick', [4 8 16 32]);
            xlabel(ax, 'MB meridians  (SBM B-scans 11 / 21 / 41 / 81)'); ylabel(ax, mlab{q});
            title(ax, sprintf('%s - %s', mlab{q}, filters(d)));
            if d == 1 && q == 1, legend(ax, 'Location', 'best'); end
        end
    end
    title(tl, 'SBM B-scan sweep vs Multimeridian MB (285 us, analogous spatial densities)', 'FontWeight', 'bold');
    MMB_SaveFigure(fig, fullfile(S, 'Metrics_vs_Bscans_SBM_vs_MB'));

    a = Agg(Agg.Filter == "noDF" & Agg.Method == "2D PG-TV", :);
    fig = figure('Color', 'w', 'Position', [40 40 800 480]);
    b_ = bar(categorical(string(a.NB) + " B / " + string(a.MB_Meridians) + " mer"), [a.AcqTime_s a.MB_AcqTime_s]);
    b_(1).DisplayName = 'SBM';  b_(2).DisplayName = 'Multimeridian MB';
    ylabel('Acquisition time [s]'); legend('Location', 'northwest'); grid on;
    title('Acquisition time at analogous spatial densities');
    MMB_SaveFigure(fig, fullfile(S, 'AcqTime_SBM_vs_MB'));
end

% ---- dispersion (frequency sweep, 21 B-scans): every period, 285 us = test mean
Disp = table();
periods = unique(arrayfun(@(r) r.Exp.Period_us, R(arrayfun(@(r) r.Exp.NB == B.TestNB, R))));
for P = periods(:)'
    for fl = filters
        for me = methods
            sel = T.Period_us == P & T.NB == B.TestNB & T.Filter == fl & T.Method == me;
            Disp = [Disp; table(P, 1e6/P, fl, me, nnz(sel), mean(T.Median_mps(sel)), ...
                std(T.Median_mps(sel)), mean(T.Ep_pct(sel)), 'VariableNames', ...
                {'Period_us', 'f_Hz', 'Filter', 'Method', 'nReps', 'c_mps', 'c_sd_mps', 'Ep_pct'})]; %#ok<AGROW>
        end
    end
end
Fit = table();
if numel(periods) >= 3
    for fl = filters
        for me = methods
            a = Disp(Disp.Filter == fl & Disp.Method == me, :);
            [p, Sfit] = polyfit(a.f_Hz/1e3, a.c_mps, 1);
            R2 = 1 - Sfit.normr^2/sum((a.c_mps - mean(a.c_mps)).^2);
            Fit = [Fit; table(fl, me, height(a), p(2), p(1), R2, 'VariableNames', ...
                {'Filter', 'Method', 'nPoints', 'c0_mps', 'dcdf_mps_per_kHz', 'R2'})]; %#ok<AGROW>
        end
    end
    writetable(Disp, fullfile(S, 'Dispersion_Points.csv'));
    writetable(Fit, fullfile(S, 'Dispersion_Fits.csv'));
    writetable(Disp, fullfile(S, 'Dispersion.xlsx'), 'Sheet', 'Points');
    writetable(Fit, fullfile(S, 'Dispersion.xlsx'), 'Sheet', 'Fits');
    MB = [];
    if exist(B.MBSummary, 'file'), MB = load(B.MBSummary, 'Disp'); end
    fig = figure('Color', 'w', 'Position', [40 40 1400 560]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
    for d = 1:2
        ax = nexttile(tl); hold(ax, 'on');
        for me = methods
            a = sortrows(Disp(Disp.Filter == filters(d) & Disp.Method == me, :), 'f_Hz');
            errorbar(ax, a.f_Hz/1e3, a.c_mps, a.c_sd_mps, '-o', 'LineWidth', 1.6, 'DisplayName', "SBM " + me);
            if ~isempty(MB)
                b = sortrows(MB.Disp(MB.Disp.Filter == filters(d) & MB.Disp.Method == me, :), 'f_Hz');
                errorbar(ax, b.f_Hz/1e3, b.c_mps, b.c_sd_mps, '--s', 'LineWidth', 1.2, 'DisplayName', "MB 8 mer " + me);
            end
        end
        yline(ax, B.Ref.Theory_cR, 'k:', 'c_R quasi-static', 'HandleVisibility', 'off');
        hold(ax, 'off'); grid(ax, 'on');
        xlabel(ax, 'Excitation frequency [kHz]'); ylabel(ax, 'Phase speed [m/s]');
        title(ax, sprintf('Dispersion - %s', filters(d))); legend(ax, 'Location', 'best');
    end
    title(tl, 'SBM frequency sweep (37 x 21, N_F = P, N_P = 4) vs Multimeridian MB', 'FontWeight', 'bold');
    MMB_SaveFigure(fig, fullfile(S, 'Dispersion_c_vs_f_SBM_vs_MB'));
    disp(Fit);
end

% ---- session needle position (sum of the coherence maps)
Cm = [];
for i = 1:numel(R)
    cm = R(i).Info.CohMap;
    if isempty(cm), continue; end
    if isfield(R(i).Info, 'Detection') && ~R(i).Info.Detection.Detected, continue; end   % noise only
    if isempty(Cm), Cm = cm; Cm.C = zeros(size(cm.C)); end
    if isequal(size(cm.C), size(Cm.C)), Cm.C = Cm.C + cm.C; end
end
Session = struct();
if ~isempty(Cm)
    [~, k] = max(Cm.C(:));  [iy, ix] = ind2sub(size(Cm.C), k);
    Session.Source_mm = [Cm.x(ix) Cm.y(iy)];
    fprintf(['Session needle position (summed coherence maps): (%.2f, %.2f) mm -> set ' ...
        'Cfg.Source.Fixed_mm = [%.2f %.2f] in the single-experiment script if the needle was not moved.\n'], ...
        Session.Source_mm, Session.Source_mm);
    fig = figure('Color', 'w', 'Position', [40 40 620 560]);
    imagesc(Cm.x, Cm.y, Cm.C/max(Cm.C(:))); axis image xy; colormap(turbo); colorbar;
    hold on; plot(Session.Source_mm(1), Session.Source_mm(2), 'w+', 'MarkerSize', 14, 'LineWidth', 2); hold off;
    xlabel('x_0 (mm)'); ylabel('y_0 (mm)'); title('Summed cylindrical-fit coherence (all acquisitions)');
    MMB_SaveFigure(fig, fullfile(S, 'Session_NeedlePosition'));
end

save(fullfile(S, 'Summary_SBM.mat'), 'T', 'Rep', 'Agg', 'Disp', 'Fit', 'Session', 'B', 'Exps');
fprintf('Summary saved (%d experiments).\n', numel(R));
end

%% ----------------------------------------------------------------------- %
function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
