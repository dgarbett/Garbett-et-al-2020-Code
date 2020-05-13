function [bgmask] = getBGMask(im_in)
% getMask returns a binary mask based on image and minCellSize
% Uses a histogram-based thresholding
% Arnold Hayer 131119
% get mask
imSmooth=imfilter(im_in,fspecial('gaussian',5,2),'symmetric');

PixMinMax=double([0 round(prctile(imSmooth(:),99))]);
IntStep=ceil((PixMinMax(2)-PixMinMax(1))/300);  
[f,xi]=ksdensity(imSmooth(:),PixMinMax(1):IntStep:PixMinMax(2));
 %figure; plot(xi,f); hold on;   
[pks,locs]=findpeaks(f);
%filter peaks greater than 0.005
log_pks=pks>0.0005; %was 0.001;
pks=pks(log_pks); locs=locs(log_pks);

%plot(xi(locs),pks,'k^','markerfacecolor',[1 0 0]);hold off;
x_bgMax=xi(locs(1)); % picks the first peak (x-value of first peak)

[~,ind]=find(f>(0.01*pks(1)),1); % returns the first value of f greater than 1% of its max
x_1pct=xi(ind); % returns the corresponding intensity value
bgWidth=x_bgMax-x_1pct; % estimates the width of the background peak
threshSeg=(x_bgMax+2*bgWidth); % multiplied by 2 because of good fg/bg separationthreshold level based on the background peak
%figure, showImagesWithLinkedAxes({imSum,imerode(imSumSmooth>threshSeg,strel('disk',1))});

bgmask=imdilate(imSmooth>threshSeg,strel('disk',2));
%figure; imagesc(bgmask);
end

