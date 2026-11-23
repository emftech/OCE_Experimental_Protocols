function [DeltaPhi] = Loupas2D_Fast(ALine_Spec,M)

    F1=ALine_Spec;

    IQ = F1(1:end,1:end-1);
    I = real(IQ); Q = imag(IQ);
    IQ1 = F1(1:end,2:end);
    I1 = real(IQ1); Q1 = imag(IQ1);
    
    Num1 = (Q.*I1)-(I.*Q1);
    Den1 = (I.*I1)+(Q.*Q1);
    
    Num2 = ((Q(1:end-1,:).*I(2:end,:))-(I(1:end-1,:).*Q(2:end,:))) +...
           ((Q1(1:end-1,:).*I1(2:end,:))-(I1(1:end-1,:).*Q1(2:end,:)));
    Den2 = ((I(1:end-1,:).*I(2:end,:))+(Q(1:end-1,:).*Q(2:end,:))) +...
           ((I1(1:end-1,:).*I1(2:end,:))+(Q1(1:end-1,:).*Q1(2:end,:)));    

    Win1 = ones(M,1);
    ConvNum1 = conv2(Win1,[1],Num1,'valid');
    ConvDen1 = conv2(Win1,[1],Den1,'valid');
    
    Win2 = ones(M-1,1);
    ConvNum2 = conv2(Win2,[1],Num2,'valid');
    ConvDen2 = conv2(Win2,[1],Den2,'valid');
    
    DeltaPhi = (atan2(ConvNum1,ConvDen1)./(1+(atan2(ConvNum2,ConvDen2)/(2*pi))));

%     figure
%     imagesc(DeltaPhi)

end