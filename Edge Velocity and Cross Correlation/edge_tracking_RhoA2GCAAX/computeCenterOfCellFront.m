function [indFrontCenter,frontMap,backMap]=computeCenterOfCellFront(protvals)
%This function takes a matrix of protrusion/retraction values (each row is
%a parametrized point on the cell edge and each column is a timepoint), and
%computes the index (row number) of the center of the cell front for each
%time point

% ----- Parameters -----
timeSmoothingParam=3;
spaceSmoothingParam=11;
frontWidth=round(0.1*size(protvals,1));
backWidth=round(0.25*size(protvals,1));

% ----- Method -----
numFrames=size(protvals,2);
numPts=size(protvals,1);

% zero pad the edges of the protrusion matrix in the time dimension to
% avoid smoothing the end of the timecourse into the beginning
protvalsZeroPadded=[zeros(size(protvals,1),timeSmoothingParam) protvals zeros(size(protvals,1),timeSmoothingParam)];
smoothingFilter=ones(timeSmoothingParam,spaceSmoothingParam); smoothingFilter=smoothingFilter/sum(smoothingFilter(:));

protSmoothed=imfilter(protvalsZeroPadded,smoothingFilter,'circular');

% remove the extra columns from both ends
protSmoothed=protSmoothed(:,(timeSmoothingParam+1):(end-timeSmoothingParam));
%imagesc(protSmoothed,[-5 5]); colormap(blue_yellow); hold on;
%imagesc(bwdist(protSmoothed<0))

% compute the center of the front as the point that is the farthest from
% the nonprotruding region of the cell
distMat=bwdist(repmat(protSmoothed<0,[3 1]));   %This is to get the right distances for points near the edge of the protrusion matrix
distMat=distMat((numPts+1):(2*numPts),:);
[~,indFrontCenter]=max(distMat,[],1);

frontMap=false(size(protvals));
backMap=false(size(protvals));
for i=1:numFrames
    indsFront=(indFrontCenter(i)-frontWidth):(indFrontCenter(i)+frontWidth);
    indsFront=1+mod(indsFront-1,numPts);
    frontMap(indsFront,i)=true;
    indBackCenter=1+mod(indFrontCenter-1+round(numPts/2),numPts);
    indsBack=(indBackCenter(i)-backWidth):(indBackCenter(i)+backWidth);
    indsBack=1+mod(indsBack-1,numPts);
    backMap(indsBack,i)=true;
end
