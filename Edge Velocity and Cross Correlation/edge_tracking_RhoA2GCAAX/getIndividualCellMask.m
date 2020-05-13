function [maskFinal cellCoors] = getIndividualCellMask(image, minCellSize)
% getMask returns a binary mask based on image and minCellSize
% Uses a histogram-based thresholding assuming there are 2 levels in the
% image
% Stores centroid and area of detected objects
% Anjali Bisaria 150824, 


logimCFP=log(1+image);
normlogimCFP=mat2gray(logimCFP);
autothreshold1=multithresh(normlogimCFP,2);
automask1=normlogimCFP>autothreshold1(2)& image < 30000;
areas=cell2mat(struct2cell(regionprops(automask1,'area')))';
mask1clean=bwareaopen(automask1,1000);
maskFinal=imfill(mask1clean, 'holes'); % Fill holes in mask
% if nargin<2
%     minCellSize=4000;
% end
% % Subtract scaled blurred image to sharpen outline
% imSmooth1=imfilter(image,fspecial('disk',1),'replicate');
% image_lpf=imfilter(imSmooth1,fspecial('disk',30),'symmetric');
% imSmooth2=imSmooth1-image_lpf*0.33; % arbitrary scaling factor
% 
% % Histogram-based threshold determination
% % find range of pixel values
% PixMinMax=double([0 round(prctile(imSmooth2(:),99))]);
% IntStep=ceil((PixMinMax(2)-PixMinMax(1))/300);  
% [f,xi]=ksdensity(imSmooth2(:),PixMinMax(1):IntStep:PixMinMax(2));
% 
% %figure; plot(xi,f); hold on;
%     
% [pks,locs]=findpeaks(f);
% %filter peaks greater than 0.005
% log_pks=pks>0.0005; 
% 
% pks=pks(log_pks); locs=locs(log_pks);
% 
% %plot(xi(locs),pks,'k^','markerfacecolor',[1 0 0]);hold off;
% x_bgMax=xi(locs(1)); % picks the first peak (x-value of first peak)
% 
% [~,ind]=find(f>(0.01*pks(1)),1); % returns the first value of f greater than 0.5% of its max
% x_1pct=xi(ind); % returns the corresponding intensity value
% bgWidth=x_bgMax-x_1pct; % estimates the width of the background peak
% threshSeg=(x_bgMax+bgWidth); % multiplied by 2 because of good fg/bg separationthreshold level based on the background peak
% %figure, showImagesWithLinkedAxes({imSum,imerode(imSumSmooth>threshSeg,strel('disk',1))});
% 
% %mask_init=imerode(imSmooth2>threshSeg,strel('disk',2));
% mask_init=imSmooth2>threshSeg;
% %mask_holes_filled=imfill(mask_init,'holes');
% % With or without filling holes
% %mask=bwareaopen(mask_holes_filled,minCellSize);
% 
% 
% % Display result
% %subplot(2,2,1); imshow(mat2gray(image));
% %subplot(2,2,2);imshow(mat2gray(mask));
% 
% %maskinv=imcomplement(mask);
% seD = strel('diamond',1);% was 2
% mask = imerode(mask_init,seD);
% mask_holes_filled=imfill(mask,'holes');
% maskFinal=bwareaopen(mask_holes_filled,minCellSize);

% BWoutline = bwperim(maskFinal);
% Segout = image;
% Segout(BWoutline) = max(image(:));
% figure, imagesc(Segout);

% get centroid coordinates
objects=regionprops(maskFinal,'BoundingBox','Area','Centroid','PixelIdxList','PixelList');
objects=objects(arrayfun(@(x) x.Area>minCellSize,objects));

cellCoors(:,1)=vect(arrayfun(@(x) x.Centroid(1),objects));
cellCoors(:,2)=vect(arrayfun(@(x) x.Centroid(2),objects));
cellCoors(:,3)=vect(arrayfun(@(x) x.Area,objects));
end




