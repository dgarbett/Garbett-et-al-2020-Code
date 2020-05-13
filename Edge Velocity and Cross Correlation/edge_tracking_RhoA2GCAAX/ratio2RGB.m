function [ imOut ] = ratio2RGB(imIn, limits)
% I=ratio2RGB(imIn,[amin amax])converts a FRET Ratio matrix (typically in the range 0 to 3)
% to an RGB image
% 140417 Arnold Hayer
imIn(isnan(imIn))=0;
if nargin<2
    im1=mat2gray(imIn);
else
    im1=mat2gray(imIn,[limits(1) limits(2)]);
    
end
im2=im2uint16(im1);
cmap=parula_black(65536);
imOut=ind2rgb(im2,cmap);
end

