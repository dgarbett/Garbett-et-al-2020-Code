function [ outlined_mask ] = DrawMaskOutline( image,mask )
% DrawMaskOutline draws a green outline around masked objects in a greyscale image
% based on a mask onto an image
% Arnold, 150521

ImageScaled=imadjust(mat2gray(image));
BWoutline = bwperim(mask);
redChan=ImageScaled; redChan(BWoutline)=0;
greenChan=ImageScaled; greenChan(BWoutline)=1;
blueChan=ImageScaled; blueChan(BWoutline)=0;
outlined_mask(:,:,1)=redChan;
outlined_mask(:,:,2)=greenChan;
outlined_mask(:,:,3)=blueChan;
end

