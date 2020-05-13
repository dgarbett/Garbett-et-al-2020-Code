function labelMask=getWindowLabelMap(cellMask,windowCoors,edgeDepthDist)
% This function creates a label matrix mapping points near the edge of the
% cell to the closest parametrized edge point (specified in windowCoors).

ringMask=cellMask-imerode(cellMask,strel('disk',edgeDepthDist));
edgeCoorMask=false(size(cellMask));
edgeCoorPixInd=windowCoors(:,1)+(windowCoors(:,2)-1)*size(cellMask,1); % 
edgeCoorMask(edgeCoorPixInd)=true;
[D,L]=bwdist(edgeCoorMask);
labelMask=L;
for i=1:size(windowCoors,1)
    labelMask(L==edgeCoorPixInd(i))=i;
end
labelMask(ringMask==0)=0;
