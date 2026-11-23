function [Frames, Axes] = SBM_TwinOCT_Reorganize(Data, Info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (implementation with Claude Code)
% Date of creation (DD/MM/YY): 17-10-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBM_TwinOCT_Reorganize
%
% Sorts SBM A-lines, in the order the ATS9360 stored them, into spatial
% frames. The .bin file layout is the MB layout:
%   A-lines per B-frame (= R) x B-frames (= C) x No. Cluster (= N_F)
% so each frame is the contiguous block of C*R A-lines that covers the
% raster body, the turn samples and the start of the flyback.
%
% Inputs:
%   Data : [Nz x Ntot] (any per-A-line quantity: complex OCT signal at
%          some depths, Doppler phase, surface displacement...) or a vector
%          of length Ntot = R*C*N_F (Info.TotalAlines).
%   Info : second output of SBM_TwinOCT_Waveform (or the saved sidecar).
%
% Outputs:
%   Frames : [N_B x N_A x N_F]           if Data is a vector
%            [Nz x N_B x N_A x N_F]      otherwise
%            Row 1 = top of the raster (Vertical Max), column 1 = Horizontal
%            Min. Serpentine reversal is already undone (spatial order).
%   Axes   : x/y axes [V] (and [mm] if Cfg.HVC/VVC are set), frame times.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isVec = isvector(Data);
if isVec, Data = reshape(Data, 1, []); end
Nper = Info.R * Info.C;
nF   = Info.Cfg.nFrames;
if size(Data,2) < Nper*nF
    error('SBM:Size', 'Data has %d A-lines, %d expected (R*C*N_F).', ...
        size(Data,2), Nper*nF);
end
Nz = size(Data,1);
D  = reshape(Data(:,1:Nper*nF), Nz, Nper, nF);
lin = Info.PixelIdx(:);                                 % column-major over (b,col)
Frames = D(:, lin, :);
Frames = reshape(Frames, Nz, Info.Cfg.nBscans, Info.Cfg.nAscans, nF);
if isVec, Frames = reshape(Frames, Info.Cfg.nBscans, Info.Cfg.nAscans, nF); end

Axes = struct();
Axes.x_V = Info.xAxis_V;  Axes.y_V = Info.yAxis_V;
if ~isempty(Info.Cfg.HVC), Axes.x_mm = Info.xAxis_V / Info.Cfg.HVC; end
if ~isempty(Info.Cfg.VVC), Axes.y_mm = Info.yAxis_V / Info.Cfg.VVC; end
Axes.t_equiv_s = (0:nF-1)' * Info.Ts_SBM;               % equivalent time (Eq. 56)
Axes.t_true_s  = (Info.FrameStart - Info.FrameStart(1)) * Info.Cfg.Ts_OCT;
end
