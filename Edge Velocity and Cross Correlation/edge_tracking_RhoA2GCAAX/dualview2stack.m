function imOut=dualview2stack(im)

[nrows ncols]=size(im);

%im1=im(:,round(0.005*ncols):round(0.4784*ncols));
im1=im(:,round(0.01*ncols):round(0.48*ncols));
%im2=im(:,round(0.4885*ncols):round(0.9770*ncols));
im2=im(:,round(0.51*ncols):round(0.98*ncols));

ncolsOut=max(size(im1,2),size(im2,2));

imOut=nan(nrows,ncolsOut,2);
imOut(:,1:size(im1,2),1)=im1;
imOut(:,1:size(im2,2),2)=im2;
