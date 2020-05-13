function [edgeCoors,edgeCoorsSmoothed]=parametrizeCellEdge(cellMask,numPoints,pointDensitySmoothingParameter,previousEdgeCoors)
%This function takes a binary cell mask and generates a matrix of
%coordinates for numPoints equally spaced points on the cell edge. If
%previousEdgeCoors is given, the output coordinates are phased and aligned
%to best match the previous coordinates.

if nargin<2
    numPoints=300;
end
if nargin<3
    pointDensitySmoothingParameter=0;
end

if pointDensitySmoothingParameter>0
    %Initially compute extraPointFactor times as many points in order to get good
    %matching of points from frame to frame
    edgeCoorsSmoothed=parametrizeCellEdge(imclose(cellMask,strel('disk',pointDensitySmoothingParameter)),numPoints,0);
    edgeCoors=mapSmoothedEdgeCoorsToCellBoundary(edgeCoorsSmoothed,cellMask);
else
    boundary=bwboundaries(cellMask);
    if length(boundary)>1
        fprintf('Warning - more than one cell detected. Using only the first one.\n');
    end
    boundary=boundary{1};
    boundaryLength=size(boundary,1)-1;
    indexCol=vect(0:(numPoints/boundaryLength):numPoints);
    edgeCoors(:,1)=interp1(indexCol,boundary(:,1),0.5:numPoints);
    edgeCoors(:,2)=interp1(indexCol,boundary(:,2),0.5:numPoints);
    edgeCoorsSmoothed=edgeCoors;
end

%  figure; hold on;
%  imagesc(cellMask);
%  scatter(edgeCoors(:,2),edgeCoors(:,1),3,'MarkerEdgeColor','g');

if nargin>3
    if size(previousEdgeCoors,1)~=numPoints
        error('Mismatch in the number of parametrized edge points');
    end
    distSum=zeros(numPoints-1,1);
    for i=1:(numPoints-1)
        temp=matrixCircularlyPermuteRows(edgeCoors,i-1);
        temp=temp(1:end,:);
        distSum(i)=sum(vect(temp-previousEdgeCoors).^2);
    end
    [~,phaseVal]=min(distSum);
    edgeCoors=matrixCircularlyPermuteRows(edgeCoors,phaseVal-1);
    edgeCoorsSmoothed=matrixCircularlyPermuteRows(edgeCoorsSmoothed,phaseVal-1);
end
