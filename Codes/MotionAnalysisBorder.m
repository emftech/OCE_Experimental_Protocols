function [loaded_phases] = MotionAnalysisBorder(Cplx_Matrix,Border,Options)


[m,n,o]=size(Cplx_Matrix);
% Options.LineOrFrame = 1; % Line = 1, Frame = 2
% Options.MotionType = 1; % Displacement = 1, particle volicity = 2;
% Options.LoupasAxialWin = 10; % 10 pixels along depth, if empty, no Loupas;
% Options.SmoothingWinPer = 0.01; % if empty, no smoothing 

%Mask = Options.Bmask;
Delta = Options.Delta;
raw_phases_Smth = [];
raw_phases_Diff = [];

 
        if Options.MotionType == 1 %%%%%%%%%%%%%%%%%%%% Displacement
            
            loaded_phases = zeros(m,o);
            if isempty(Options.SmoothingWinPer) % No smoothing
                
                %%%%%%%%%%%%%%%%%
                for ii = 1:m
                    if isnan(Border(ii))
                        loaded_phases(ii,1:end) = NaN(1,o);
                    else
                        raw_phases = -1*angle(squeeze(mean(Cplx_Matrix(ii,Border(ii):Border(ii)+Delta,:),2)));
                        loaded_phases(ii,1:end) = raw_phases';
                    end
                    ii
                end
                %%%%%%%%%%%%%%%%
            
            else % Smoothing
                
                %%%%%%%%%%%%%%%%%
                for ii = 1:m
                    if isnan(Border(ii))
                        loaded_phases(ii,1:end) = NaN(1,o);
                    else
                        raw_phases = -1*angle(squeeze(mean(Cplx_Matrix(ii,Border(ii):Border(ii)+Delta,:),2)));
                        raw_phases_Smth = smooth(raw_phases,Options.SmoothingWinPer,'lowess')';

                        loaded_phases(ii,1:end) = raw_phases_Smth';
                    end
                    ii
                end
                %%%%%%%%%%%%%%%%
            end
            
        elseif Options.MotionType == 2 %%%%%%%%%%%% particle volicity
            if isempty(Options.LoupasAxialWin) % No Loupas
                
                loaded_phases = zeros(m,o-1);
                if isempty(Options.SmoothingWinPer) % No smoothing

                    %%%%%%%%%%%%%%%%%
                    for ii = 1:m
                        if isnan(Border(ii))
                            loaded_phases(ii,1:end) = NaN(1,o-1);
                        else
                            raw_phases = -1*angle(squeeze(mean(Cplx_Matrix(ii,Border(ii):Border(ii)+Delta,:),2)));
                            raw_phases_Diff = diff(unwrap(raw_phases,[],1),1,1);
                            loaded_phases(ii,1:end) = raw_phases_Diff';
                        end
                        ii
                    end
                    %%%%%%%%%%%%%%%%                    
                        
                else % Smoothing
                    
                    %%%%%%%%%%%%%%%%%
                    for ii = 1:m
                        if isnan(Border(ii))
                            loaded_phases(ii,1:end) = NaN(1,o-1);
                        else
                            %idx = find(Mask(Border(ii):Border(ii)+Delta,ii)==1);
                            %raw_phases = -1*angle(squeeze(mean(Cplx_Matrix(ii,Border(ii)+idx-1,:),2)));
                                                        
                            raw_phases = -1*angle(squeeze(mean(Cplx_Matrix(ii,Border(ii):Border(ii)+Delta,:),2)));
                            raw_phases_Diff = diff(unwrap(raw_phases,[],1),1,1);
                            raw_phases_Smth = smooth(raw_phases_Diff,Options.SmoothingWinPer,'lowess')';

                            loaded_phases(ii,1:end) = raw_phases_Smth';
                        end
                        %ii
                    end
                    %%%%%%%%%%%%%%%%
                
                end
            
            else % Loupas
                
                loaded_phases = zeros(m,o-1);
                if isempty(Options.SmoothingWinPer) % No smoothing
                    
                    %%%%%%%%%%%%%%%%%
                    for ii = 1:m
                        raw_phases_Diff = Loupas2D_Fast(squeeze(Cplx_Matrix(ii,:,:)),Options.LoupasAxialWin);
                        if isnan(Border(ii)) | (Border(ii)-round(Options.LoupasAxialWin/2))<=0
                            loaded_phases(ii,1:end) = NaN(1,o-1);
                        else
                            loaded_phases(ii,1:end) = raw_phases_Diff(Border(ii)-round(Options.LoupasAxialWin/2),:);
                        end
                        ii
                    end
                    %%%%%%%%%%%%%%%%                    

                else % Smoothing

                    %%%%%%%%%%%%%%%%%
                    for ii = 1:m
                        raw_phases_Diff = Loupas2D_Fast(squeeze(Cplx_Matrix(ii,:,:)),Options.LoupasAxialWin);
                        if isnan(Border(ii)) | (Border(ii)-round(Options.LoupasAxialWin/2))<=0
                            loaded_phases(ii,1:end) = NaN(1,o-1);
                        else
                            %raw_phases_Diff_line = raw_phases_Diff(Border(ii)-round(Options.LoupasAxialWin/2),:);
                            raw_phases_Diff_line = mean(raw_phases_Diff(Border(ii)-round(Options.LoupasAxialWin/2):...
                                                                        Border(ii)-round(Options.LoupasAxialWin/2)+Delta,:),1);
                            raw_phases_Smth = smooth(raw_phases_Diff_line,Options.SmoothingWinPer,'lowess')';
                            loaded_phases(ii,1:end) = raw_phases_Smth;
                        end
                      
                        ii
                    end
                    %%%%%%%%%%%%%%%%                    
                                       
                end
                
            end           
        
        end     
           
        
           
        
end




