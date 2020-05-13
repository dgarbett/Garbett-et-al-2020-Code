clear;clc;close all; 
root = 'K:\data\150626_NG';
rawdir=[root,filesep,'CC-Rac1_1min'];
bgdir=[root,filesep,'background'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[rawdir,filesep,'data_150711'];
position='Pos2';
%%
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
bgrawFRET=double(imread([bgdir,filesep,'AVG_bgEPI-FRET.tif']));
bgFRET1=dualview2stack(bgrawFRET);
bgFRET2=dualviewAlignFromFittedSurface(bgFRET1,pX,pY,binning);
CFPbg=bgFRET2(:,:,2);
FRETbg=bgFRET2(:,:,1);
folder=[rawdir,filesep,position];
%folder=[rawdir,filesep,position,filesep,'Pos0'];
filenames=getFilenames(folder);
filenamesFRET=filenames(boolRegExp(filenames,'FRET'));

imRatio_raw={};maskFinal={};cellCoors={};
for frameNum=11%:length(filenamesFRET)
    disp(num2str(frameNum));
    imrawFRET=double(imread([folder,filesep,filenamesFRET{frameNum}]));
    imFRET1=dualview2stack(imrawFRET);
    imFRET2=dualviewAlignFromFittedSurface(imFRET1,pX,pY,binning);
     
    imCFP_raw=imFRET2(:,:,2);
    imFRET_raw=imFRET2(:,:,1);
   
    %%%%%% Align CFP/FRET images
    
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imCFP_raw+imFRET_raw);
    imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
    imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
    %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMask((imCFP_raw),2000);
    maskFinal{frameNum}=mask;
    cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Detrmine ratio
    imFRETbg(~mask)=nan;
    imFRET=ndnanfilter(imFRETbg,fspecial('disk',1),'replicate');
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',1),'replicate');
    imRatioTemp=imFRET./imCFP;
    imRatioTemp(~mask)=nan;
    imRatio_raw{frameNum}=imRatioTemp;
    %%%%%% Determine scaling for representation
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),1),1),round(prctile(imRatioTemp(:),99),1)];
    end
    colorRange=[0.9 1.1];
    %%%%%% Generate and write files for raw ratio and outlined objects
    tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);
   
    %imwrite(imFRETOutline,[datadir,filesep,position,'_FRET.tif'],'WriteMode','append','Compression','none');
    %imwrite(tempRATIO,[datadir,filesep,position,'_RATIO_raw_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
    %imwrite(imRFPOutline,[datadir,filesep,position,'_RFP.tif'],'WriteMode','append','Compression','none');
end


    
%% Display alignment
figure;
subplot(1,2,1);showImagesMergeChannels(imFRET1(:,:,1),imFRET1(:,:,2));
subplot(1,2,2);showImagesMergeChannels(imCFP_raw,imFRET_raw);

%% Display segmentation
figure;
testimage=[imFRETOutline{frameNum} tempRATIO];
imshow(testimage);

%% background mask
figure; imagesc(bgmask);