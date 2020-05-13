function imOut=dualview2stack_AH(im)

[nrows ncols]=size(im);

imOut(:,:,1)=im(:,481:960);

imOut(:,:,2)=im(:,1:480);


