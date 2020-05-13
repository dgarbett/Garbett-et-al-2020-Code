function [maskFinal,cellCoors]=fretGetCellMasks(minSize, imCFP, imFRET, mask_raw)
%fretGetCellMasks takes as input an image and a set of bounding rectangles
%(output from fretGetInitialObjectsAndBGmask) and computes a more precise
%mask for all cells detected except those touching the image boundary

if nargin<3 || isempty(minSize)
    minSize=2000;
end


maskFinal=imfilter(mask_raw,fspecial('disk',2),'symmetric');
%figure; imagesc(maskFinal);

objects=regionprops(maskFinal,'BoundingBox','Area','Centroid','PixelIdxList','PixelList');
objects=objects(arrayfun(@(x) x.Area>minSize,objects));

cellCoors(:,1)=vect(arrayfun(@(x) x.Centroid(1),objects));
cellCoors(:,2)=vect(arrayfun(@(x) x.Centroid(2),objects));
cellCoors(:,3)=vect(arrayfun(@(x) x.Area,objects));
for i=1:length(objects)
    cellCoors(i,4:7)=objects(i).BoundingBox;
    if nargin>3
        pixInd=objects(i).PixelIdxList;
        cellCoors(i,10)=sum(imFRET(pixInd))/sum(imCFP(pixInd));
        cellCoors(i,11)=median(imFRET(pixInd)./imCFP(pixInd));
        
        %CFP and FRET intensities
        cellCoors(i,12)=mean(imCFP(pixInd));
        cellCoors(i,13)=mean(imFRET(pixInd));
        
        %Smoothed values for estimating polarization direction
        cfpSmooth=nan(size(imFRET)); cfpSmooth(pixInd)=imCFP(pixInd);
        cfpSmooth=ndnanfilter(cfpSmooth,fspecial('gaussian',7,2),7);
        fretSmooth=nan(size(imFRET)); fretSmooth(pixInd)=imFRET(pixInd);
        fretSmooth=ndnanfilter(fretSmooth,fspecial('gaussian',7,2),7);
        ratioSmooth=fretSmooth./cfpSmooth;
        
        %Weighted centroid:
        cellCoors(i,8)=sum(objects(i).PixelList(:,1).*ratioSmooth(pixInd))/sum(ratioSmooth(pixInd));
        cellCoors(i,9)=sum(objects(i).PixelList(:,2).*ratioSmooth(pixInd))/sum(ratioSmooth(pixInd));
        
        %Measures of polarization
        polarizationVector=[cellCoors(i,8)-cellCoors(i,1) cellCoors(i,9)-cellCoors(i,2)];
        polarizationVector=polarizationVector/norm(polarizationVector);
        relCoord=objects(i).PixelList - repmat(cellCoors(i,1:2),[length(pixInd) 1]);
        projectedCoord=relCoord*polarizationVector(:);
        cellCoors(i,14)=myNanCorrcoef(projectedCoord,ratioSmooth(pixInd));
        ptile75=prctile(projectedCoord,75);
        ptile25=prctile(projectedCoord,25);
        pixRear=pixInd(projectedCoord<=ptile25);
        pixFront=pixInd(projectedCoord>=ptile75);
        cellCoors(i,15)=sum(imFRET(pixRear))/sum(imCFP(pixRear));
        cellCoors(i,16)=sum(imFRET(pixFront))/sum(imCFP(pixFront));
    end
end

%------------------------------------------------------------------------------
function maskFinal=processSingleObjectRectangle(maskFinal,objRect,imForSeg0,minSize)

width=fretGetObjectRectanglePaddingSize;

cx1=floor(max(1,objRect(1)));
cy1=floor(max(1,objRect(2)));
cx2=ceil(min(size(imForSeg0,2),sum(objRect([1 3]))));
cy2=ceil(min(size(imForSeg0,1),sum(objRect([2 4]))));

%Smooth and sharpen the image before segmenting
imForSeg=imForSeg0(cy1:cy2,cx1:cx2);
imForSeg=imfilter(imForSeg,fspecial('gaussian',5,2)); %Try to change// MASK!!!
imSumNorm=mat2gray(imForSeg, [prctile(imForSeg(:),5) prctile(imForSeg(:),99.5)]);
imSumNorm=imsharpenSean(imSumNorm,4,30); %Try enhancing the edges//MASK!!! changed from 50 to 30 130827

%Segment the image
threshSeg=graythresh(imSumNorm(imSumNorm>0 & imSumNorm<1));
maskThisObj=imSumNorm>threshSeg;
maskThisObj=bwareaopen(maskThisObj,25,4); %showImagesWithLinkedAxes({imForSeg0(cy1:cy2,cx1:cx2),imSumNorm,imSumNorm>threshSeg,maskThisObj});
maskThisObj=imdilate(maskThisObj,strel('disk',1));
maskThisObj=imclose(maskThisObj,strel('disk',2));
%maskThisObj=imclearborder(maskThisObj);
maskThisObj=imfill(maskThisObj,'holes'); %plotCDFs({imSumNorm(:),[threshSeg]});
maskThisObj=bwareaopen(maskThisObj,minSize);

objFound=regionprops(maskThisObj,'BoundingBox');

if length(objFound)>1
    maskFinal(cy1:cy2,cx1:cx2)=maskFinal(cy1:cy2,cx1:cx2) | maskThisObj;
else
    maskFinal(cy1:cy2,cx1:cx2)=maskFinal(cy1:cy2,cx1:cx2) | maskThisObj;
end