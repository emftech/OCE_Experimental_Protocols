function [Frame_tmp_r,X,Y] = Elastogram(ParamsElasto,OCT_system,OCE_system,Frames2,Bmode_Mask,Bmode_IntLog,Border,Xaxis,Zaxis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log 16/08/2026:
% Edits to Fernando's original Elastogram function in order for the function to work: 
% 1. fixed interp2 in Speed2D to SpeedEquivInterp conversion (~line 32)
% 2. Added "omitnan" on meridian mean speed estimation (~lines 145 158, and 159)
% 3. Changed "gcf" to "f" in figure creation and saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
Params = ParamsElasto;
Ts = Params.Ts;

freq = linspace(0,1,Params.FFT_Nx)/Ts;
Nf = find(freq>Params.Freq);
Nf = Nf(1);

Params.FFT_FreqSample = Nf;
DepthCut = Params.ZcutRange;


%% Wave speed 2D - multi-meridian

[Speed2D, K, Xaxis_E, Zaxis_E] = SpeedEstimation_PhaseDeriv(Frames2(1:end,DepthCut(1):DepthCut(2),:),Xaxis(1:end),Zaxis(1:diff(DepthCut)+1),Params);
Speed2D_Equiv = abs(Speed2D.Equiv);
Speed2D_Equiv(Speed2D_Equiv>20) = NaN;
Speed2D_Equiv(Speed2D_Equiv<0) = NaN;
Speed2D_Equiv_Filt = nanmedfilt2(Speed2D_Equiv,[3 3]);

%Speed2D_Equiv_Filt = medfilt2(Speed2D_Equiv,[3 3],'Symmetric');

figure
imagesc(Xaxis_E,Zaxis_E,abs(Speed2D_Equiv_Filt))
caxis([0 20])
colormap(jet)

% SpeedEquivInterp = interp2(Xaxis_E,Zaxis_E,Speed2D_Equiv_Filt,Xaxis,Zaxis(1:diff(DepthCut)+1),'linear',0.0005);
[Xq,Zq] = meshgrid(Xaxis, Zaxis(1:diff(DepthCut)+1));
SpeedEquivInterp = interp2( ...
    Xaxis_E, ...
    Zaxis_E, ...
    Speed2D_Equiv_Filt, ...
    Xq, ...
    Zq, ...
    'linear', ...
    0.0005);


% figure
% imagesc(Xaxis,Zaxis,SpeedEquivInterp.*Bmode_Mask)
% caxis([0 5])

f = figure;
imagesc(Xaxis,Zaxis(1:diff(DepthCut)+1),SpeedEquivInterp.*Bmode_Mask(DepthCut(1):DepthCut(2),:))
ylabel('z-axis (mm)');
xlabel('x-axis (mm)');
title(['Speed Map - ',num2str(Params.Freq),' Hz'])
%axis equal
axis([Xaxis_E(1) Xaxis_E(end) 0 Zaxis(diff(DepthCut)+1)])
set(gca,'FontSize',14)
colorbar
colormap(jet)
caxis([0 20])
saveas(f,['Speed_Map_Alone_',num2str(Params.Freq),'Hz_OLD.fig']);
saveas(f,['Speed_Map_Alone_',num2str(Params.Freq),'Hz_OLD.tif']);


I = uint8(255*mat2gray(Bmode_IntLog(DepthCut(1):DepthCut(2),:).*...
                       Bmode_Mask(DepthCut(1):DepthCut(2),:), ...
                       [OCE_system.i_thresh_low OCE_system.i_thresh_high]));
bg = ind2rgb(I,gray(255));
bgImg = double(bg);
alphaFactor = 0.4;
bgImgAlpha_in = (1 - alphaFactor) .* bgImg.*Bmode_Mask(DepthCut(1):DepthCut(2),:);
bgImgAlpha_out = bgImg.*(1-Bmode_Mask(DepthCut(1):DepthCut(2),:));
bgImgAlpha = bgImgAlpha_in + bgImgAlpha_out;

II = uint8(255*mat2gray(SpeedEquivInterp.*Bmode_Mask(DepthCut(1):DepthCut(2),:), [0 20]));
im = ind2rgb(II,hot(255));

fgImg = double(im);
fgImgAlpha = alphaFactor .* fgImg.*(Bmode_Mask(DepthCut(1):DepthCut(2),:));
fusedImg = fgImgAlpha + bgImgAlpha;

f = figure;
imagesc(Xaxis,Zaxis(1:diff(DepthCut)+1),fusedImg)
%imagesc(Xaxis1*1e3,Yaxis1*1e3,squeeze(SpaceTime(:,:,j)))
ylabel('z-axis (mm)');
xlabel('x-axis (mm)');
%axis equal
axis([Xaxis_E(1) Xaxis_E(end) 0 Zaxis(diff(DepthCut)+1)])
title(['Speed Map - ',num2str(Params.Freq),' Hz'])
set(gca,'FontSize',14)
colorbar
colormap(hot)
caxis([0 20])
drawnow
saveas(f,['Elastogram_2D_Hot_',num2str(Params.Freq),'Hz_OLD.fig']);
saveas(f,['Elastogram_2D_Hot_',num2str(Params.Freq),'Hz_OLD.tif']);


I = uint8(255*mat2gray(Bmode_IntLog(DepthCut(1):DepthCut(2),:).*Bmode_Mask(DepthCut(1):DepthCut(2),:), [OCE_system.i_thresh_low OCE_system.i_thresh_high]));
bg = ind2rgb(I,gray(255));
bgImg = double(bg);
alphaFactor = 0.4;
bgImgAlpha_in = (1 - alphaFactor) .* bgImg.*Bmode_Mask(DepthCut(1):DepthCut(2),:);
bgImgAlpha_out = bgImg.*(1-Bmode_Mask(DepthCut(1):DepthCut(2),:));
bgImgAlpha = bgImgAlpha_in + bgImgAlpha_out;

II = uint8(255*mat2gray(SpeedEquivInterp.*Bmode_Mask(DepthCut(1):DepthCut(2),:), [0 20]));
im = ind2rgb(II,jet(255));

fgImg = double(im);
fgImgAlpha = alphaFactor .* fgImg.*(Bmode_Mask(DepthCut(1):DepthCut(2),:));
fusedImg = fgImgAlpha + bgImgAlpha;


f = figure;
imagesc(Xaxis,Zaxis(1:diff(DepthCut)+1),fusedImg)
%imagesc(Xaxis1*1e3,Yaxis1*1e3,squeeze(SpaceTime(:,:,j)))
ylabel('z-axis (mm)');
xlabel('x-axis (mm)');
%axis equal
axis([Xaxis_E(1) Xaxis_E(end) 0 Zaxis(diff(DepthCut)+1)]) 
title(['Speed Map - ',num2str(Params.Freq),' Hz'])
set(gca,'FontSize',14)
colorbar
colormap(jet)
caxis([0 20])
drawnow
saveas(f,['Elastogram_2D_Jet1_',num2str(Params.Freq),'Hz_OLD.fig']);
saveas(f,['Elastogram_2D_Jet1_',num2str(Params.Freq),'Hz_OLD.tif']);

%% En-face wave speed map

NumBscanSingle = OCT_system.NumBscanSingle;
DeltaSpace = Params.DeltaSpace;
Center_idx = NumBscanSingle/2;

for i = 1:OCT_system.NumAngles
    
    Speed_tmp = SpeedEquivInterp(:,NumBscanSingle*(i-1)+1:NumBscanSingle*(i))';
    
    for j = 1: NumBscanSingle 
        pointer_ini = Border.Idx_Up(NumBscanSingle*(i-1)+j)-DepthCut(1)+1;
        if pointer_ini<1 | pointer_ini > diff(DepthCut)+1-DeltaSpace-10
            pointer_ini = round(diff(DepthCut)/2);
        end
        % Speed_Line(j) = median(Speed_tmp(j,pointer_ini:pointer_ini+DeltaSpace),2); 
        Speed_Line(j) = median(Speed_tmp(j,pointer_ini:pointer_ini+DeltaSpace),2,"omitnan"); % added omitmissing to avoid big patches
    end
    
%     figure
%     imagesc(Speed_tmp')
%     colormap(jet)
%     caxis([0 20])
    
    SpeedMer(i*2-1,:) = Speed_Line(NumBscanSingle/2:-1:1);
    SpeedMer(i*2,:) = Speed_Line(NumBscanSingle/2+1:NumBscanSingle);
    
    % SpeedMerAve(i*2-1,1) = mean(Speed_Line(1,Center_idx-DeltaSpace-4:Center_idx-4));
    % SpeedMerAve(i*2,1) = mean(Speed_Line(1,Center_idx+4:Center_idx+4+DeltaSpace));
    SpeedMerAve(i*2-1,1) = mean(Speed_Line(1,Center_idx-DeltaSpace-4:Center_idx-4),"omitnan"); % added omitmissing to avoid big patches
    SpeedMerAve(i*2,1) = mean(Speed_Line(1,Center_idx+4:Center_idx+4+DeltaSpace),"omitnan"); % added omitmissing to avoid big patches


end

AnglesAll = [0:180/OCT_system.NumAngles:180/OCT_system.NumAngles*(2*OCT_system.NumAngles-1)]';
[idx_angles] = IDXSelection(OCT_system.NumAngles);

SpeedMer_New = SpeedMer(idx_angles,:,:);
SpeedMerAve_New = SpeedMerAve(idx_angles,:,:);

theta = repmat(AnglesAll/180*pi,1,NumBscanSingle/2);
rho = repmat(linspace(0,OCT_system.x_scan_width/2,NumBscanSingle/2),length(AnglesAll),1);
[x,y] = pol2cart(theta,rho);

[X, Y] = meshgrid(linspace(-OCT_system.x_scan_width/2,OCT_system.x_scan_width/2,NumBscanSingle));

% figure
% scatter(x,y)
    
rhoSpeed = SpeedMer_New(:,:);

F = scatteredInterpolant(x(:),y(:),rhoSpeed(:),'natural','none');
Frame_tmp_r = F(X,Y);
Frame_tmp_r = medfilt2(Frame_tmp_r,[5 5],'Symmetric');

f = figure;
imagesc(X(1,:),Y(:,1),Frame_tmp_r)
%imagesc(Xaxis1*1e3,Yaxis1*1e3,squeeze(SpaceTime(:,:,j)))
ylabel('y-axis (mm)');
xlabel('x-axis (mm)');
axis equal
axis([-OCT_system.x_scan_width/2 OCT_system.x_scan_width/2,...
      -OCT_system.x_scan_width/2 OCT_system.x_scan_width/2]);
colormap(jet)
colorbar
caxis([0 20]) 
set(gca,'FontSize',14)
drawnow
saveas(f,['Elastogram_EnFace_Jet1_',num2str(Params.Freq),'Hz_OLD.fig']);
saveas(f,['Elastogram_EnFace_Jet1_',num2str(Params.Freq),'Hz_OLD.tif']);

% uiopen('FigModelPolarPlot.fig',1)
% gcf = figure(1);
f = figure;
f.Position = [481 339 837 351];
subplot(1,2,1)
polarplot([AnglesAll/180*pi; 0],[SpeedMerAve_New; SpeedMerAve_New(1)],'-bs','linewidth',1,'MarkerSize',8)
title(['Phase speed (m/s) @ ',num2str(Params.Freq),'Hz.'])
set(gca,'FontSize',14);
subplot(1,2,2)
boxplot(SpeedMerAve_New,'Whisker',1);
[Mean,Error,Range] = grpstats (SpeedMerAve_New,[],{'mean','std','range'});
title(['Mean ',num2str(round(Mean,2)),' m/s, STD ',num2str(round(Error,3)),' m/s.'])
grid on
saveas(f,'PolarPhaseSpeedAvg_OLD.fig');
saveas(f,'PolarPhaseSpeedAvg_OLD.tif');

%% Saving variables

save('Elastograms_OLD.mat','Speed2D','Xaxis_E','Zaxis_E','Speed2D_Equiv_Filt',...
                      'SpeedEquivInterp','SpeedMer_New','SpeedMerAve_New',...
                      'theta','rho','X','Y','rhoSpeed','Frame_tmp_r');

                  
end
