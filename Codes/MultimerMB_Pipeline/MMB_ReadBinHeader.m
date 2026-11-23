function [ScanInfo, headerBytes, Check] = MMB_ReadBinHeader(binFile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Estéfano Machado Ferrer
% Date of creation: 07-10-2026
%
% MMB_ReadBinHeader
%
% Reads the text header of the LabVIEW .bin OCE files (same fields as
% oct2_readRawData) and checks that the file holds the complete raw volume
%   2 bytes x samples x A-lines x B-frames x meridians.
%
% Outputs:
%   ScanInfo    : header fields (samples_in_Aline, Alines_in_Bframe,
%                 Bframes_in_3Dscan, No_3Dscans, scan lengths, ...)
%   headerBytes : offset of the raw uint16 data
%   Check       : ExpectedBytes, ActualBytes, IsComplete, and
%                 CompleteMeridians (meridians fully present on disk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(binFile);
if fid < 0
    error('Cannot open %s', binFile);
end
L = fread(fid, 1, 'uint16');
fseek(fid, 4, 'bof');
H = fread(fid, L, '*char')';
fclose(fid);

headerBytes = 4 + L;

st = strfind(H, ': ') + 2;
en = strfind(H, char(13)) - 3;
v = cell(1, numel(st)-1);
for c = 2:numel(st)
    v{c-1} = H(st(c):en(c));
end

ScanInfo = struct();
ScanInfo.acquisition_datetime = v{1};
ScanInfo.sample_rate_MHz = str2double(v{3});
ScanInfo.samples_in_Aline = str2double(v{4});
ScanInfo.type = v{7};
ScanInfo.Alines_in_Bframe = str2double(v{8});
ScanInfo.Flyback_Aline_No = str2double(v{9});
ScanInfo.Bframe_trigger_Aline_delay = str2double(v{10});
ScanInfo.Bframes_in_3Dscan = str2double(v{11});
ScanInfo.No_3Dscans = str2double(v{12});
ScanInfo.Hor_scan_length_mm = str2double(v{13});
ScanInfo.Ver_scan_length_mm = str2double(v{14});

d = dir(binFile);
bytesPerMer = 2*ScanInfo.samples_in_Aline*ScanInfo.Alines_in_Bframe*ScanInfo.Bframes_in_3Dscan;

Check = struct();
Check.ExpectedBytes = headerBytes + bytesPerMer*ScanInfo.No_3Dscans;
Check.ActualBytes = d.bytes;
Check.IsComplete = d.bytes >= Check.ExpectedBytes;
Check.CompleteMeridians = floor((d.bytes - headerBytes)/bytesPerMer);

end
