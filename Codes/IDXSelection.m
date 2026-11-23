function [idx] = IDXSelection(NumAngles)

Jumps1 = repmat([3 1],1,NumAngles);
Jumps1 = Jumps1(1:NumAngles-1);
Jumps2 = repmat([1 3],1,NumAngles);
Jumps2 = Jumps2(1:NumAngles-1);

idx1 = cumsum([1, Jumps1]);
idx2 = cumsum([2, Jumps2]);

idx = [idx1 idx2];

end