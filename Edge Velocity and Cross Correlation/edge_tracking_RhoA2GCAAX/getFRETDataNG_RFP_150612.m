function getFRETDataNG(position)
% image processing for FRET data analysis 
% position='Pos0';
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

%% Set up
root = 'K:\data\150603_NG';
rawdir=[root,filesep,'RhoA2G_CDH5-mRuby3_timelapse_1min\RhoA2G_CDH5-mRuby3_1min'];
bgdir=[root,filesep,'background'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[root,filesep,'RhoA2G_CDH5-mRuby3_timelapse_1min\data_150611_1'];
if ~exist(datadir)
    mkdir(datadir);
end

%% %%%% Call background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
bgrawFRET=double(imread([bgdir,filesep,'AVG_bgEPI-FRET.tif']));
bgFRET1=dualview2stack(bgrawFRET);
bgFRET2=dualviewAlignFromFittedSurface(bgFRET1,pX,pY,binning);
CFPbg=bgFRET2(:,:,2);
FRETbg=bgFRET2(:,:,1);

bgrawRFP=double(imread([bgdir,filesep,'AVG_bgEPI-RFP.tif']));
bgRFP1=dualview2stack(bgrawRFP);
bgRFP2=dualviewAlignFromFittedSurface(bgRFP1,pX,pY,binning);
RFPbg=bgRFP2(:,:,1);

%%%%%% Loop through frames
folder=[rawdir,filesep,position];
filenames=getFilenames(folder);
filenamesFRET=filenames(boolRegExp(filenames,'FRET'));
filenamesRFP=filenames(boolRegExp(filenames,'RFP'));

imRatio_raw={};maskFinal={};cellCoors={};
for frameNum=1:length(filenamesFRET)
    disp(num2str(frameNum));
    imrawFRET=double(imread([folder,filesep,filenamesFRET{frameNum}]));
    imFRET1=dualview2stack(imrawFRET);
    imFRET2=dualviewAlignFromFittedSurface(imFRET1,pX,pY,binning);
    
    imrawRFP=double(imread([folder,filesep,filenamesRFP{frameNum}]));
    imRFP1=dualview2stack(imrawRFP);
    imRFP2=dualviewAlignFromFittedSurface(imRFP1,pX,pY,binning);
    
    imCFP_raw=imFRET2(:,:,2);
    imFRET_raw=imFRET2(:,:,1);
    imRFP_raw=imRFP2(:,:,1);
    %%%%%% Align CFP/FRET images
    
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imCFP_raw);
    imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
    imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
    %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMask((imCFP_raw),2000);
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
       colorRange = [round(prctile(imRatioTemp(:),3),1),round(prctile(imRatioTemp(:),97),1)];
    end
    %%%%%% Generate and write files for raw ratio and outlined objects
    %tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);
    imRFPOutline{frameNum}=DrawMaskOutlineRuby(imRFP_raw,mask);
    %imwrite(imFRETOutline,[datadir,filesep,position,'_FRET.tif'],'WriteMode','append','Compression','none');
%     imwrite(tempRATIO,[datadir,filesep,position,'_RATIO_raw_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
    %imwrite(imRFPOutline,[datadir,filesep,position,'_RFP.tif'],'WriteMode','append','Compression','none');
end
%% %%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay
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
        colorRange = [round(prctile(frame_1(:),5),1),round(prctile(frame_1(:),97),1)];
%        colorRange = [3.5, 5.5];
   end
   tempRATIO=imRatio_raw{frameNum}./corr(frameNum);
   imRatio{frameNum}=tempRATIO;
   tempRATIOforstack=ratio2RGB(tempRATIO,colorRange);%Cdc42
   stitched=[imFRETOutline{frameNum} imRFPOutline{frameNum} tempRATIOforstack];
   imwrite(stitched,[datadir,filesep,position,'_stitched_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
   disp(num2str(frameNum));
end       
%%%%%% Save data
save([datadir,filesep,position,'_RatioData.mat'],'maskFinal','cellCoors','imRatio');
end

