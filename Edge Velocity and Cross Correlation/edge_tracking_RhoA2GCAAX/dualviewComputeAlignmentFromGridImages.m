function [pX,pY,dxMat,dyMat] = dualviewComputeAlignmentFromGridImages(imGrid)

% Prepare images for aligning
ims=dualview2stack(imGrid);
im1g=double(imfilter(ims(:,:,1),fspecial('gaussian',5,2),'symmetric')); 
im2g=double(imfilter(ims(:,:,2),fspecial('gaussian',5,2),'symmetric'));
%im1g=imgradient(im1g);
%im2g=imgradient(im2g);

% Compute empirical local alignments
num=30; %was 15
step=floor(size(im1g,1)/num);
blockWidth=60;
ctrsX=round(step/2):step:size(im1g,2);
ctrsY=round(step/2):step:size(im1g,1);
dxMat=zeros(length(ctrsY),length(ctrsX));
dyMat=dxMat;
for i=1:length(ctrsY)
    disp(num2str(i));
    for j=1:length(ctrsX)
        y0=max(ctrsY(i)-blockWidth,1); y1=min(ctrsY(i)+blockWidth,size(im1g,1));
        x0=max(ctrsX(j)-blockWidth,1); x1=min(ctrsX(j)+blockWidth,size(im1g,2));
        [dxMat(i,j), dyMat(i,j)]=getOffsetBetweenImages(im1g(y0:y1,x0:x1),im2g(y0:y1,x0:x1));
    end
end

% Fit a smooth function to the empirical data - I will use a second order
% surface
xMat=dualviewGetXMatGridFor2ndOrderSurface(ctrsX,ctrsY);

pX=robustfit(xMat,dxMat(:));
pY=robustfit(xMat,dyMat(:));
