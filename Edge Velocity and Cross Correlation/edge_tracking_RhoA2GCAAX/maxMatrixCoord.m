function [i1,i2] = maxMatrixCoord(mat)
% Finds the coordinates of the maximum value in a matrix

[maxV,maxInd]=max(mat);
[~,i2]=max(maxV);
i1=maxInd(i2);
