function [ScanInfo, headerBytes, Check] = SBM_ReadBinHeader(binFile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer (with Claude Code)
% Date of creation: 19-10-2026
%
% SBM_ReadBinHeader
%
% Reads the text header of a LabVIEW .bin file (same fields as
% MMB_ReadBinHeader / oct2_readRawData, parsed by key instead of position)
% and checks that the file holds the complete SBM acquisition:
%   2 bytes x samples x (A-lines per B-frame = R) x (B-frames = C) x (3D = N_F)
%
% Outputs:
%   ScanInfo    : samples_in_Aline, Alines_in_Bframe (R), Bframes_in_3Dscan (C),
%                 No_3Dscans (N_F), Bframe_trigger_Aline_delay, scan lengths,
%                 Flyback_Aline_No (N_turn), acquisition_datetime, Raw
%   headerBytes : offset of the raw uint16 data
%   Check       : ExpectedBytes, ActualBytes, IsComplete, CompleteFrames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(binFile, 'r');
if fid < 0
    error('SBM:File', 'Cannot open %s', binFile);
end
L = fread(fid, 1, 'uint16');
fseek(fid, 4, 'bof');
H = fread(fid, L, '*char')';
fclose(fid);
headerBytes = 4 + L;

tok = regexp(H, '([^:;\r\n]+):\s*([^;\r\n]*)', 'tokens');
key = strtrim(cellfun(@(t) t{1}, tok, 'UniformOutput', false));
val = cellfun(@(t) strtrim(t{2}), tok, 'UniformOutput', false);

ScanInfo = struct();
ScanInfo.samples_in_Aline           = getNum('spectral samples per A-line');
ScanInfo.sample_rate_MHz            = getNum('Sample Rate');
ScanInfo.Alines_in_Bframe           = getNum('A-lines per B-frame');
ScanInfo.Flyback_Aline_No           = getNum('Flyback duration');
ScanInfo.Bframe_trigger_Aline_delay = getNum('B-frame trigger delay');
ScanInfo.Bframes_in_3Dscan          = getNum('B-frames in 3D scan');
ScanInfo.No_3Dscans                 = getNum('3D repeated scans');
ScanInfo.Hor_scan_length_mm         = getNum('Horizontal scan length');
ScanInfo.Ver_scan_length_mm         = getNum('Vertical scan length');
ScanInfo.Max_Alines_per_buffer      = getNum('Max A-lines per buffer');
ScanInfo.acquisition_datetime       = getStr('Date/Time scan acquisition start');
ScanInfo.type                       = getStr('Scan type');
ScanInfo.Raw = H;
try
    ScanInfo.Datetime = datetime(ScanInfo.acquisition_datetime, ...
        'InputFormat', 'M/d/yyyy h:mm:ss a', 'Locale', 'en_US');
catch
    ScanInfo.Datetime = NaT;
end

d = dir(binFile);
bytesPerFrame = 2*ScanInfo.samples_in_Aline*ScanInfo.Alines_in_Bframe*ScanInfo.Bframes_in_3Dscan;
Check = struct();
Check.ExpectedBytes  = headerBytes + bytesPerFrame*ScanInfo.No_3Dscans;
Check.ActualBytes    = d.bytes;
Check.IsComplete     = d.bytes >= Check.ExpectedBytes;
Check.CompleteFrames = floor((d.bytes - headerBytes)/bytesPerFrame);

    function v = getNum(k)
        i = find(contains(key, k), 1);
        if isempty(i), v = NaN; else, v = str2double(val{i}); end
    end
    function s = getStr(k)
        i = find(contains(key, k), 1);
        if isempty(i), s = ''; else, s = val{i}; end
    end
end
