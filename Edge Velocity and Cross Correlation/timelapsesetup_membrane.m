function [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_membrane(mem_temp,nucr,debrisarea,badframes,RGBwrite,RGBdir)
%%% determine median cell size for blur detection %%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=zeros(3,1); numcells=zeros(3,1);
for i=1:3
    raw1=mem_temp(:,:,i);
    %nuc_mask=blobdetector(raw1,nucr,blobthreshold,debrisarea);
    mem_mask=threshmask(raw1,3);
    mem_mask=markershed(mem_mask,round(nucr*2/3));
    mem_mask=bwareaopen(mem_mask,debrisarea);
   
    [~,numcells(i)]=bwlabel(mem_mask);
    mem_area(i)=median(cell2mat(struct2cell(regionprops(mem_mask,'Area'))));
end
dims=size(mem_mask);
height=dims(1); width=dims(2);
blurthreshhigh=1.10*nanmedian(mem_area);
%%% determine first good frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstgoodindex=1;

badframes(1)=0;
blurthreshhigh=1.5*mem_area(firstgoodindex);
blurthreshlow=0.5*mem_area(firstgoodindex);
numthresh=0.5*numcells(firstgoodindex); %was 0.5
end

