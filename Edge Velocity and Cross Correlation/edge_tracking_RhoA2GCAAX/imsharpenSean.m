function imSharp = imsharpenSean(im,magnitude,gaussianRadius,threshold)
%IMSHARPENSEAN implements image sharpening
%   Detailed explanation goes here

imBlurred=imfilter(im,fspecial('gaussian',max(3,round(2.5*gaussianRadius)),gaussianRadius),'symmetric');
sharpMask=im-imBlurred; %imagesc(sharpMask); colorbar;
imSharp=im+magnitude*sharpMask;
%imSharp=im-magnitude*imBlurred;
%showImagesWithLinkedAxes({im,imSharp});

end

