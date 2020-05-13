function edgeCoors=mapSmoothedEdgeCoorsToCellBoundary(ecSmooth,cellMask,weightDistToNormal)
% This function maps a smoothed parametrized cell edge boundary (or any
% parametrized edge boundary) to nearby points on the boundary of a
% cellMask. By default, it minimizes a weighted average of the distance to
% the point on the parametrized boundary and a distance to normal line
% intersecting that point.

if nargin<3
    weightDistToNormal=2;
end

boundary=bwboundaries(cellMask);
if length(boundary)>1
    fprintf('Warning - more than one cell detected. Using only the first one.\n');
end
boundary=boundary{1};

normV=computeNormalVectorsFromParametrizedCellEdge(ecSmooth,round(size(ecSmooth,1)/12));
%nnInd = annsearch(boundary', ecSmooth', 1);
nnInd=zeros(size(ecSmooth,1),1);
for i=1:length(nnInd)
    distances=weightDistToNormal*distanceFromPointsToLine(boundary,[ecSmooth(i,:) normV(i,:)]) + sqrt(sum((boundary-repmat(ecSmooth(i,:),[size(boundary,1) 1])).^2,2));
    [~,temp]=min(distances);
    nnInd(i)=temp(1);
end

edgeCoors=boundary(nnInd,:);