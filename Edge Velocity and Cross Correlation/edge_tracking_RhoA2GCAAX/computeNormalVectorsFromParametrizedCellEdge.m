function normalVectors=computeNormalVectorsFromParametrizedCellEdge(edgeCoors,numPtsToUse)
% This function takes a set of parametrized coordinates for a cell edge and
% computes a normal vector for each point

if nargin<2
    numPtsToUse=7;
end

numBoundaryPts=size(edgeCoors,1);

ecWrapped=[edgeCoors((end-numPtsToUse):end,:); edgeCoors; edgeCoors(1:numPtsToUse,:)];
diffmat=diff(ecWrapped);
diffmat=0.5*(diffmat(1:(end-1),:) + diffmat(2:end,:));  %Use the symmetric diferences (n+1) - (n-1)
diffmatSmoothed=[smooth(diffmat(:,1),numPtsToUse) smooth(diffmat(:,2),numPtsToUse)];

gradientMat=diffmatSmoothed(numPtsToUse:(end-numPtsToUse),:);
normalVectors=[-1*gradientMat(:,2) gradientMat(:,1)]./repmat(sqrt(gradientMat(:,1).^2 + gradientMat(:,2).^2),[1 2]);
