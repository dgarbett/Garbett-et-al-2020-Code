function getFRETData( row,col,site )
% image processing for FRET data analysis 
%
% Input: CFP/FRET images, background images, and alignment parameters for
% CFP/FRET images. 
%
% Generates first a mask for background subtraction, subtracts background,
% and calculates a refined mask before taking a ratio between background-
% subtracted CFP and FRET images. Bleaching correction assumes constant 
% mean FRET values per frame across the entire time-lapse series.
%
% Output: 
%
%   maskFinal   Cell array with binary mask for each frame
%
%   cellCoors   Cell array with centroid x,y coordinates of detected
%               objects and object size in pixels
%
%   imRatioraw  Cell array with masked raw ratio values, before bleaching
%               correction
%
%   imRatio     Cell array with bleaching-corrected masked ratio values
%
%   imFRETOUtline_row_col.tif 
%               RGB tif stack of FRET channel images with outlined detected
%               objects
%
%   imRatio_row_col_colorRangeLow_colorRangeHigh.tif
%               RGB tif stack with false-colored (jet-black) scaled ratio
%               images, bleaching corrected.
%
% Used non built-in subfunctions: 
%
%   dualviewAlignFromFittedSurface
%   getBGMask
%   subBG
%   getCellMask
%   DrawMaskOutline
%   ratio2RGB
%   vect
%
% Arnold, 24 May 2015

%%%%%% Set up
root = 'K:\data\150424_IX\RhoFRET-CDH5-CHY-2min';
rawdir=[root,filesep,'rawdata'];
bgdir='K:\data\150424_IX\RhoFRET-CDH5-CHY-2min\background-bin1';
load([root,filesep,'alignment parameters pX pY.mat'],'pX','pY');
datadir=[root,filesep,'data_150521'];
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%%%%%% Call background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
CFPbg_raw=double(imread([bgdir,filesep,'AVG_bgCFP.tif']));
FRETbg_raw=double(imread([bgdir,filesep,'AVG_bgFRET.tif']));
bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);
CFP_files=dir([rawdir,filesep,shot,'_CFP_*']);
%%%%%% Loop through frames
imRatio_raw={};maskFinal={};cellCoors={};
for frameNum=1:length(CFP_files)
    disp([num2str(row),'_',num2str(col),'_',num2str(site),'_',(num2str(frameNum))]);
    imCFP_raw=double(imread([rawdir,filesep,shot,'_CFP_',num2str(frameNum),'.tif']));
    imFRET_raw=double(imread([rawdir,filesep,shot,'_FRET_',num2str(frameNum),'.tif']));
    %%%%%% Align CFP/FRET images
    imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
    imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
    imCFP_raw=imaligned(:,:,1);
    imFRET_raw=imaligned(:,:,2);
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imCFP_raw+imFRET_raw);
    imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
    imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
    %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMask(imFRET_raw,2000);
    maskFinal{frameNum}=mask;
    cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Detrmine ratio
    imFRETbg(~mask)=nan;
    imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
    imRatioTemp=imFRET./imCFP;
    imRatioTemp(~mask)=nan;
    imRatio_raw{frameNum}=imRatioTemp;
    %%%%%% Determine scaling for representation
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),2),1),round(prctile(imRatioTemp(:),98),1)];
    end
    %%%%%% Generate and write files for raw ratio and outlined objects
    tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imFRETOutline=DrawMaskOutline(imFRET_raw,mask);
    imwrite(imFRETOutline,[datadir,filesep,shot,'_FRET.tif'],'WriteMode','append','Compression','none');
    imwrite(tempRATIO,[datadir,filesep,shot,'_RATIO_raw_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
end
%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay
bleach_1=nanmean(vect(imRatio_raw{1}));
for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
end
bleach_corr=bleach_raw/bleach_1;
timepts=1:length(imRatio_raw);
[fitvals]=polyfit(timepts,bleach_corr,1);
corr=polyval(fitvals,timepts);
%plot(timepts,polyval(fitvals,timepts));
%%%%%% Correct imRatio_raw for bleaching
imRatio={};
for frameNum=1:length(imRatio_raw)
   if frameNum==1
       frame_1=imRatio_raw{1};
       colorRange = [round(prctile(frame_1(:),5),1),round(prctile(frame_1(:),95),1)];
   end
   imRatio{frameNum}=imRatio_raw{frameNum}./corr(frameNum);
   tempRATIO=ratio2RGB(imRatio{frameNum},colorRange);%Cdc42
   imwrite(tempRATIO,[datadir,filesep,shot,'_RATIO_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
   %disp(num2str(frameNum));
end       
%%%%%% Save data
save([datadir,filesep,shot,'_RatioData.mat'],'maskFinal', 'cellCoors','imRatio_raw', 'imRatio');
end
