function Results = OCE_PG_ProcessMB(SpaceTime, fs, f_exc, scanLength_mm, nPosPerMer, nMer, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
%
% OCE_PG_ProcessMB
%
% Runs the four phase-speed estimation strategies of the thesis on one
% EXPERIMENTAL Multimeridian MB-mode acquisition:
%
%   1D PG    : 1D phase gradient along each meridian, then en face
%              interpolation of the speed estimates.
%   1D PG-TV : 1D phase gradient, en face interpolation of |k|, 2D TV on |k|.
%   2D PG    : en face interpolation of every sparse frame (Delaunay),
%              then 2D phase gradient (local plane fit).
%   2D PG-TV : 2D phase gradient followed by TV on |k|.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   SpaceTime     : [nMer*nPosPerMer x Nt] motion traces, one row per
%                   lateral position in acquisition order (e.g.
%                   loaded_phases_Border from MotionAnalysisBorder)
%   fs            : A-line rate [Hz] (OCT_system.a_scan_rate*1e3)
%   f_exc         : excitation frequency [Hz]
%   scanLength_mm : meridian length [mm] (ScanInfo.Ver_scan_length_mm)
%   nPosPerMer    : positions per meridian (ScanInfo.Bframes_in_3Dscan)
%   nMer          : number of meridians (ScanInfo.No_3Dscans)
%
% Optional name-value parameters:
%   Geometry : 'ReverseAlternate' (true), 'AngleOffset_deg' (0),
%              'FlipRadial' (false), 'Center_mm' ([0 0])
%   Shared   : 'KernelSize_px' (11), 'R2Threshold' (0.7),
%              'SpeedLimits' ([0.5 30]), 'AmpThreshold' (0.05),
%              'ExcludeCenterRadius_mm' (0), 'GridSize' (nPosPerMer),
%              'PhasorMethod' ('dtft'), 'TemporalHann' (true)
%   TV       : 'MuTV' ([] -> MuTVRel*median|k|), 'MuTVRel' (0.1),
%              'TVIterations' (500), 'TVTau' (0.25), 'TVSigma' (0.25),
%              'TVTheta' (1), 'TVBoundary' ('neumann')
%   1D       : 'InterpMethod1D' ('natural')
%   2D       : 'InterpMethod2D' ('linear'), 'Smoothing2D' (true),
%              'Sigma_px' (0.8)
%   Control  : 'Run1D' (true), 'Run2D' (true), 'KeepMovie' (false),
%              'Verbose' (true)
%
% Output struct Results:
%   Geom, Meridian_Frames, PG1D, PG2D, Interp, (Movie), Params, and a
%   convenience struct Results.Maps with the four masked elastograms:
%       Maps.PG1D, Maps.PG1D_TV, Maps.PG2D, Maps.PG2D_TV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.FunctionName = 'OCE_PG_ProcessMB';

addParameter(p, 'ReverseAlternate', true);
addParameter(p, 'AngleOffset_deg', 0);
addParameter(p, 'FlipRadial', false);
addParameter(p, 'Center_mm', [0 0]);

addParameter(p, 'KernelSize_px', 11);
addParameter(p, 'R2Threshold', 0.70);
addParameter(p, 'SpeedLimits', [0.5 30]);
addParameter(p, 'AmpThreshold', 0.05);
addParameter(p, 'ExcludeCenterRadius_mm', 0);
addParameter(p, 'GridSize', nPosPerMer);
addParameter(p, 'PhasorMethod', 'dtft');
addParameter(p, 'TemporalHann', true);

addParameter(p, 'MuTV', []);
addParameter(p, 'MuTVRel', 0.1);
addParameter(p, 'TVIterations', 500);
addParameter(p, 'TVTau', 0.25);
addParameter(p, 'TVSigma', 0.25);
addParameter(p, 'TVTheta', 1.0);
addParameter(p, 'TVBoundary', 'neumann');

addParameter(p, 'InterpMethod1D', 'natural');
addParameter(p, 'InterpMethod2D', 'linear');
addParameter(p, 'Smoothing2D', true);
addParameter(p, 'Sigma_px', 0.8);

addParameter(p, 'Run1D', true);
addParameter(p, 'Run2D', true);
addParameter(p, 'KeepMovie', false);
addParameter(p, 'Verbose', true);

parse(p, varargin{:});
opts = p.Results;

%% ===================== Geometry and meridian frames ======================

Geom = OCE_MB_MeridianGeometry(nPosPerMer, nMer, scanLength_mm, ...
    'ReverseAlternate', opts.ReverseAlternate, ...
    'AngleOffset_deg', opts.AngleOffset_deg, ...
    'FlipRadial', opts.FlipRadial, ...
    'Center_mm', opts.Center_mm);

Meridian_Frames = OCE_MB_ReshapeToMeridians(SpaceTime, nPosPerMer, nMer);

% Sanity check of the excitation frequency: c = 2*pi*f/|k| is directly
% proportional to the frequency used by the estimator.
[f_dom, FreqInfo] = OCE_EstimateDominantFrequency(SpaceTime, fs);
if abs(f_dom - f_exc) > FreqInfo.Resolution_Hz/2
    warning(['Dominant motion frequency (%.1f Hz) differs from f_exc (%.1f Hz) by more ' ...
             'than half the spectral resolution (%.1f Hz). Check f_exc and the A-line rate.'], ...
             f_dom, f_exc, FreqInfo.Resolution_Hz);
end

sharedPG = { ...
    'KernelSize_px', opts.KernelSize_px, ...
    'R2Threshold', opts.R2Threshold, ...
    'SpeedLimits', opts.SpeedLimits, ...
    'AmpThreshold', opts.AmpThreshold, ...
    'ExcludeCenterRadius_mm', opts.ExcludeCenterRadius_mm, ...
    'PhasorMethod', opts.PhasorMethod, ...
    'TemporalHann', opts.TemporalHann, ...
    'MuTV', opts.MuTV, ...
    'MuTVRel', opts.MuTVRel, ...
    'TVIterations', opts.TVIterations, ...
    'TVTau', opts.TVTau, ...
    'TVSigma', opts.TVSigma, ...
    'TVTheta', opts.TVTheta, ...
    'TVBoundary', opts.TVBoundary, ...
    'Verbose', opts.Verbose};

Results = struct();
Results.Geom = Geom;
Results.Meridian_Frames = Meridian_Frames;
Results.fs = fs;
Results.f_exc = f_exc;
Results.f_dominant = f_dom;
Results.Params = opts;
Results.Maps = struct();

%% ===================== 1D PG and 1D PG-TV ================================

if opts.Run1D
    Results.PG1D = OCE_PG1D_AlongMeridians(Meridian_Frames, Geom, fs, f_exc, ...
        sharedPG{:}, ...
        'GridSize', opts.GridSize, ...
        'InterpMethod', opts.InterpMethod1D);

    Results.x_mm = Results.PG1D.x_mm;
    Results.y_mm = Results.PG1D.y_mm;
    Results.Maps.PG1D    = Results.PG1D.Elastogram_1DPG_Masked;
    Results.Maps.PG1D_TV = Results.PG1D.Elastogram_1DPG_TV_Masked;
end

%% ===================== 2D PG and 2D PG-TV ================================

if opts.Run2D
    [Movie, Interp] = OCE_MB_InterpolateEnFaceFrames(Meridian_Frames, Geom, ...
        'GridSize', opts.GridSize, ...
        'InterpMethod', opts.InterpMethod2D, ...
        'Smoothing', opts.Smoothing2D, ...
        'Sigma_px', opts.Sigma_px, ...
        'Verbose', opts.Verbose);

    Results.PG2D = OCE_PG2D(Movie, Interp.x_mm, Interp.y_mm, fs, f_exc, ...
        sharedPG{:}, ...
        'SupportMask', Interp.HullMask, ...
        'Center_mm', opts.Center_mm);

    Results.Interp = Interp;
    Results.x_mm = Interp.x_mm;
    Results.y_mm = Interp.y_mm;
    Results.Maps.PG2D    = Results.PG2D.Elastogram_2DPG_Masked;
    Results.Maps.PG2D_TV = Results.PG2D.Elastogram_2DPG_TV_Masked;

    if opts.KeepMovie
        Results.Movie = Movie;
    end
end

%% ===================== Summary ===========================================

if opts.Verbose
    fprintf('\nOCE_PG_ProcessMB summary (masked maps, m/s):\n');
    names = fieldnames(Results.Maps);
    for i = 1:numel(names)
        v = Results.Maps.(names{i});
        v = v(isfinite(v));
        fprintf('  %-8s median = %6.3f, IQR = [%6.3f %6.3f], valid px = %d\n', ...
            names{i}, median(v), prctile(v,25), prctile(v,75), numel(v));
    end
end

end
