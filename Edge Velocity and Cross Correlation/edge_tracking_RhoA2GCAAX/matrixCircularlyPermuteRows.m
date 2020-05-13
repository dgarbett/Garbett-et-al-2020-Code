function mat=matrixCircularlyPermuteRows(mat,shift)
%This function circularly permutes the rows of a matrix

shift=mod(shift,size(mat,1));
if shift<0
    mat=[mat((1-shift):end,:); mat(1:(-1*shift),:)];
end
if shift>0
    mat=[mat((end-shift+1):end,:); mat(1:(end-shift),:)];
end
