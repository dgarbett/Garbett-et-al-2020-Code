function [imOut,bg] = imageSubtractBackgroundWithObjectMaskingAndPadding(im,blockSize,padSize,objMask)
%imageSubtractBackgroundWithObjectMasking does image background subtraction
%after an initial pass of separating foreground from background elements.
%It computes block by block with the given blocksize, but also includes a
%surrounding region (determined by padsize) to get more pixel data.

if nargin<2
    blockSize=80;
end
if nargin<3
    padSize=60;
end

if nargin<4   %Automatically do an initial segmentation of background and foreground
    imNorm=mat2gray(im, [prctile(im(:),5) prctile(im(:),99.5)]);
    threshInitialSeg=graythresh(imNorm(imNorm>0 & imNorm<1));
    maskBGinitial=imNorm<threshInitialSeg;
    maskFGinitial=imNorm>threshInitialSeg;
    maskBGinitial=imdilate(maskBGinitial,strel('disk',5));
    maskBGinitial=imerode(maskBGinitial,strel('disk',20));
    %maskBGinitial=~imdilate(maskFGinitial,strel('disk',5));
    imMaskedInitial=im;
    imMaskedInitial(~maskBGinitial)=0;
else          %Use the object mask passed in
    imMaskedInitial=im;
    imMaskedInitial(~objMask)=0;
end

%Do the background subtraction
[~,bg]=imageSubtractBackgroundWithPadding(imMaskedInitial,blockSize,padSize);  %zeros are ignored in the background computation
imOut=im-bg;

end

