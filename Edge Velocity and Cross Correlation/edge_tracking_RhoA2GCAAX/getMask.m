function [mask] = getMask(image,  minCellSize )
% getMask returns a binary mask based on image and minCellSize
% Uses a histogram-based thresholding
% Arnold Hayer 131119

if nargin<2
    minCellSize=4000;
end
%imSmooth=imfilter(image,fspecial('gaussian',3,2),'symmetric');


imSmooth=imfilter(image,fspecial('disk',1),'replicate');
%imagesc(imSmooth);


[f,xi]=ksdensity(vect(imSmooth),-50:4:10000);

%figure; plot(xi,f); hold on;
    
[pks,locs]=findpeaks(f);
%filter peaks greater than 0.005
log_pks=pks>0.0002; 

pks=pks(log_pks); locs=locs(log_pks);

%plot(xi(locs),pks,'k^','markerfacecolor',[1 0 0]);hold off;
x_bgMax=xi(locs(1)); % picks the first peak (x-value of first peak)

[~,ind]=find(f>(0.005*pks(1)),1); % returns the first value of f greater than 1% of its max
x_1pct=xi(ind); % returns the corresponding intensity value
bgWidth=x_bgMax-x_1pct; % estimates the width of the background peak
threshSeg=(x_bgMax+bgWidth); % multiplied by 2 because of good fg/bg separationthreshold level based on the background peak
%figure, showImagesWithLinkedAxes({imSum,imerode(imSumSmooth>threshSeg,strel('disk',1))});

mask_init=imerode(imSmooth>threshSeg,strel('disk',2));
mask_holes_filled=imfill(mask_init,'holes');
% With or without filling holes
%mask=bwareaopen(mask_holes_filled,minCellSize);
mask=bwareaopen(mask_init,minCellSize);

% Display result
%subplot(2,2,1); imshow(mat2gray(image));
%subplot(2,2,2);imshow(mat2gray(mask));

maskinv=imcomplement(mask);
% subplot(2,2,3); 
% rgb(:,:,1)=0.8*(mat2gray(maskinv));
% rgb(:,:,2)=mat2gray(log(image));
% rgb(:,:,3)=zeros(size(mask));
% imshow(rgb);
end




