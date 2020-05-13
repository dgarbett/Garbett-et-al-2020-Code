function [dx,dy]=getOffsetBetweenImages(im1,im2,offsetRange,scaleFactor)
% Aligns to images and computes the offset in x and y between them. Some of
% the strategies are borrowed from Feng-Chiao's CalcJitter functions

if nargin<3
    offsetRange=20; %was 9
end
if nargin<4
    scaleFactor=3;  %Default is to use 1/3 pixel accuarcy
end

% Original unprocessed images
im10=double(im1);
im20=double(im2);

%Smooth both images a little bit to minimize pixel noise effects. Then
%compute the gradient of each image - this should emphasize the edges in
%the image that we want to align.
if scaleFactor>1
    im1=zeros(size(im10)*scaleFactor);
    im2=zeros(size(im20)*scaleFactor);
    for i=1:3
        for j=1:3
            im1(i:scaleFactor:end,j:scaleFactor:end)=im10;
            im2(i:scaleFactor:end,j:scaleFactor:end)=im20;
        end
    end
end
% im1=double(imfilter(im1,fspecial('gaussian',7,2.5),'symmetric')); 
% im2=double(imfilter(im2,fspecial('gaussian',7,2.5),'symmetric'));
im1=double(imfilter(im1,fspecial('gaussian',9,5),'symmetric')); 
im2=double(imfilter(im2,fspecial('gaussian',9,5),'symmetric'));
im1=imgradient(im1);
im2=imgradient(im2);

%Normalize the gradient images such that |im1 - im2| is a good similarity
%metric and to zero out a lot of the background pixels.
im1=im1-nanmedian(im1(:)); 
im2=im2-nanmedian(im2(:)); 
if (nanmean(im1(:))>0) && (nanmean(im2(:))>0)
    im1=max(im1,0)/nanmean(im1(:));
    im2=max(im2,0)/nanmean(im2(:));
end

%% preparation
ymatrix=repmat((-1*offsetRange:offsetRange)',[1 2*offsetRange+1]);
xmatrix=repmat((-1*offsetRange:offsetRange),[2*offsetRange+1 1]);

%% calculate distance scores for a grid of offset values
dmat=zeros(2*offsetRange+1,2*offsetRange+1);
for i=1:(2*offsetRange+1)
    for j=1:(2*offsetRange+1)

        dy=ymatrix(i,j);dx=xmatrix(i,j);

        im1cropped=im1(max(1,1-dy):min(end,end-dy),max(1,1-dx):min(end,end-dx));  %crop the appropriate side of the image, depending on whether dx and dy are positive or negative
        im2cropped=im2(max(1,1+dy):min(end,end+dy),max(1,1+dx):min(end,end+dx));
%        ImageD=abs(im1cropped-im2cropped);
        ImageD=(im1cropped-im2cropped).*(im1cropped-im2cropped);
        
        dmat(i,j)=nanmean(ImageD(:));  % The mean absolute difference between the images will be used as a distance metric

    end
end

%% find the offsets with the minimum distance score
[i,j]=maxMatrixCoord(-1*dmat);
minV=dmat(i,j);
dx=xmatrix(i,j)/scaleFactor;
dy=ymatrix(i,j)/scaleFactor; %[xmatrix(i,j) ymatrix(i,j)]
%showImagesMergeChannels(shiftMatrix(im1,3*dx,3*dy),im2);

%% interpolate the distance score matrix to estimate the offset between the images with subpixel accuracy
% % If desired, this section can be used to do subpixel registration
% dmatInterp=imresize(dmat,[501 501],'bicubic');
% xmatrixInterp=imresize(xmatrix,[501 501],'bicubic');
% ymatrixInterp=imresize(ymatrix,[501 501],'bicubic');
% [iInterp,jInterp]=maxMatrixCoord(-1*dmatInterp);
% dx=xmatrixInterp(iInterp,jInterp);
% dy=ymatrixInterp(iInterp,jInterp);

