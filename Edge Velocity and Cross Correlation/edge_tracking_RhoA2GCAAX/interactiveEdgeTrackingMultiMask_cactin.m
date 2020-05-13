function interactiveEdgeTrackingMultiMask(position,bgdir,rawdir,datadir,jitter) %#ok<INUSD>
%%%%meant for fret images with the beam splitter
%run bg_images_alignment_script before hand... try to build in before hand
%% Set up
filesep= '/';
%load([bgdir,filesep,'alignment parameters pX pY.mat']);
jitter=0;
jitterx=25;
jittery=60;
binning=1;
position = '';

Cy5_stack=[];FRET_stack=[];
for i = 1:40;

    shot=['Pos7b_CAAX_g003'];
    Cy5_temp=imread([rawdir,filesep,shot,'_Cy5_',sprintf('%02d',1),'.tif']);
    YFP_temp=imread([rawdir,filesep,shot,'_YFP_',sprintf('%02d',1),'.tif']);
    Cy5_stack(:,:,i)=Cy5_temp;
    YFP_stack(:,:,i)=YFP_temp;

end
Cy5_bg=uint16(mean(Cy5_stack,3));
YFP_bg=uint16(mean(YFP_stack,3));

imwrite(Cy5_bg,[bgdir,filesep,'AVG_rawdata_Cy5.tif'],'TIFF','Compression','None');
imwrite(YFP_bg,[bgdir,filesep,'AVG_rawdata_YFP.tif'],'TIFF','Compression','None');

Cy5_bg = double(Cy5_bg);
YFP_bg = double(YFP_bg);


%open images and perform the dual view alignment
folder=[rawdir,filesep,position];
filenames=getFilenames(folder);
filenamesYFP=filenames(boolRegExp(filenames,'YFP'));
filenamesCy5=filenames(boolRegExp(filenames,'Cy5'));
imRatio_raw={};maskFinal={};cellCoors={};

%load([bgdir,filesep,'alignment parameters pX pY.mat']);


%%for each position(k) open up the first image and last image
imrawRFP=double(imread([folder,filesep,filenamesRFP{1}]));
%imRFP1=dualview2stack(imrawRFP);
%imRFP2=dualviewAlignFromFittedSurface(imRFP1,pX,pY,binning);
%imRFP_raw=imRFP2(:,:,1);

imrawRFP=double(imread([folder,filesep,filenamesRFP{end}]));
%imRFP1=dualview2stack(imrawRFP);
%imRFP2=dualviewAlignFromFittedSurface(imRFP1,pX,pY,binning);
%imRFP_raw_end=imRFP2(:,:,1);
clear imRFP1 imRFP2;

figure(1)
h=imagesc(log(imrawRFP));
axis image
title('First Image')
colormap(gray);
%%
%user input the number of cells to track
prompt = {'Enter Number of Cells:'};
dlg_title = 'Input';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
numCells = str2num(answer{1});

%makes sure the user thinks there are cells.
if (numCells < 1)
    return;
end;

%populate an array with the positions of each ROI
cellRegions = cell(numCells,1);
figure(1)

for c = 1:numCells
    % ask user for another mask
    area = imfreehand(gca);
    setColor(area,'red');
    cellRegions{c} = createMask(area);
end

close all;
%go through each area and perform analysis:

cellMask = cellRegions;
cellCoors = {};
maskFinal = {};
%%
currentMap = parula;
newMap = [0 0 0; currentMap];
colormap(newMap)

for frameNum=1:length(filenamesFRET)
    
    %load in images
    imrawYFP=double(imread([folder,filesep,filenamesYFP{frameNum}]));   
    imrawCy5=double(imread([folder,filesep,filenamesCy5{frameNum}]));
 

    % jitter correction using RFP image
    if jitter==1;
        if frameNum>1
            [alignpara ~]=dftregistration(fft2(imCy5_raw0),fft2(imrawCy5),100);
            xx=alignpara(4);yy=alignpara(3);
            x=x+xx;y=y+yy;
        else
            x=0;y=0;    %store jitters data
        end
        imCy5_raw0=imrawCy5;
        disp([x,y])
        
        imrawYFP=CropJitter(imrawYFP,jitterx,jitterx,jittery,jittery,x,y);
        imrawCy5=CropJitter(imrawCy5,jitterx,jitterx,jittery,jittery,x,y);
        
        for c=1:numCells;
            cellMask{c} = CropJitter(cellMask{c},jitterx,jitterx,jittery,jittery,x,y);
        end;
        %CFPbg=CropJitter(CFPbg_raw,jitterx,jitterx,jittery,jittery,x,y);
        %FRETbg=CropJitter(FRETbg_raw,jitterx,jitterx,jittery,jittery,x,y);
        %RFPbg=CropJitter(RFPbg_raw,jitterx,jitterx,jittery,jittery,x,y);
    else
        %CFPbg=CFPbg_raw;
        %FRETbg=FRETbg_raw;
        %RFPbg=RFPbg_raw;
    end
    
    
    
    maskAll = cell(numCells,1);
    imRatioAll = cell(numCells,1);
    cellCoorsAll = [];
    maskCombined = zeros(size(imrawYFP));
    
    bgmask=getBGMask(imrawYFP+imrawCy5);
    %imCy5bg=subBG(imrawCy5,bgmask,Cy5_bg);
    %imYFPbg=subBG(imrawYFP,bgmask,YFP_bg);
    
    imCy5 = ndnanfilter(imrawCy5,fspecial('disk',1),'replicate');
    imYFP = ndnanfilter(imrawYFP,fspecial('disk',1),'replicate');    
    %imYFP = imYFP - repmat(nanmin(imYFP(:)),size(imYFP)); 
    %imCy5 = imCy5 - repmat(nanmin(imCy5(:)),size(imCy5)); 
    imRatioTemp = imYFP./imCy5;
    
    %iterate through the cells to get data for each frame
    for c = 1:numCells;
        % First mask using RFP channel
        [mask cellCoorsTemp]=getIndividualCellMask((imrawCy5.*cellMask{c}),2000);
        maskAll{c,1}  =  mask;
        maskRing = imdilate(mask, strel('disk',3));
        bgCy5Px = regionprops(maskRing, imCy5,'PixelValues');
        bgCy5Av = nanmean(bgCy5Px.PixelValues(:));
        bgYFPpx = regionprops(maskRing, imYFP,'PixelValues');
        bgYFPAv = nanmean(bgYFPpx.PixelValues(:));
        %bg subtraction???
        cellCoorsAll = [cellCoorsAll; cellCoorsTemp];
        maskCombined = maskCombined |  mask;
    end %through cells

    imRatioTemp(~maskCombined) = nan;
    
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),3)*10)/10,round(prctile(imRatioTemp(:),97)*10)/10];
    end
    %tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imFRETOutline{frameNum}=DrawMaskOutline(imrawYFP,maskCombined);
    imRFPOutline{frameNum}=DrawMaskOutlineRuby(imrawCy5,maskCombined);
    
    maskFinal{frameNum}=maskAll;
    cellCoors{frameNum}=cellCoorsAll;
    imRatio_raw{frameNum}=imRatioTemp;
    
    h=figure('visible','off');
    im = imagesc(imRatioTemp);
    colormap(newMap);
    colorbar;
    caxisLow = 0;
    caxisHigh = 2;
    caxis([caxisLow caxisHigh]);
    set(gca,'XTick',[]) % Remove the ticks in the x axis!
    set(gca,'YTick',[]) % Remove the ticks in the y axis
 
    saveas(gcf,[bgdir,filesep,shot,'_Ratio_',num2str(caxisLow),'_',num2str(caxisHigh),'_',sprintf('%02d',frameNum)],'png')
    
   %imwrite(uint16(imRatioTemp),[bgdir,filesep,shot,'_Ratio',sprintf('%02d',frameNum),'.tif'],'TIFF','Compression','None');

    
end%through frames

for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
end
save([datadir,filesep,position,'_RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw')%,'imFRETOutline','imRFPOutline');
save([datadir,filesep,position,'_Bleach_raw.mat'],'bleach_raw');

end

