function [Border] = FindSurface(Bmode,OCE_system,OCT_system)

IniDepth = 10;

PeakThresMult = OCE_system.PeakThresMult;
PeakThresWinSize = OCE_system.PeakThresWinSize; 
axial_pixel_res = OCT_system.axial_pixel_res;

Bmode1=20*log10(medfilt2(Bmode,[3 3],'symmetric'));


for ii=1:OCE_system.NumLateralPos

    Profile = Bmode1(IniDepth:end,ii);
    %Thres = mean(Profile(1:PeakThresWinSize))+std(Profile(1:PeakThresWinSize))*PeakThresMult;
    % MEAN - FERNANDO
    % Thres = mean(Profile(1:PeakThresWinSize))+PeakThresMult;
    % HARMMEAN - ALONSO
    Thres = harmmean(Profile(1:PeakThresWinSize))+PeakThresMult;

    %Thres = PeakThresMult;
    
    [pks, locs] = findpeaks(Profile,'MinPeakHeight',Thres,'NPeaks',1);  

    if isempty(locs)
        Border.Idx(ii) = NaN;
        Border.DepthPos(ii) = NaN;
        Border.Inten(ii) = NaN;
    else
        Border.Idx(ii) = locs(1)+IniDepth;
        Border.DepthPos(ii) = (locs(1)-1+IniDepth)*axial_pixel_res*1e-3; % in mm
        Border.Inten(ii) = 20*log10(Bmode(locs(1)+IniDepth,ii));
    end 
    %ii    
end

end
