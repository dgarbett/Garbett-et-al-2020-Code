function getActMyoData(position,rawdir,datadir,jitter)
% image processing for FRET data analysis 
% clear; clc
% 
% root='K:\data\150205_IXM\FT-CHY-P2A-TQ-MYL9-CDH5-CIT-2min-si29';
% position='2_6_1';
% rawdir=[root,filesep,'rawdata'];
% datadir=root;
% jitter=0;


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


%%%%%% Loop through frames
folder=rawdir;
%folder=[rawdir,filesep,position,filesep,'Pos0'];
filenames=getFilenames(folder);
filenamesRFP=filenames(boolRegExp(filenames,[position,'_Texas']));
filenamesCFP=filenames(boolRegExp(filenames,[position,'_CFP']));

imRatio_raw={};maskFinal={};cellCoors={};
for frameNum=1:length(filenamesRFP)
    
    disp(num2str(frameNum));
    imCFPraw=double(imread([folder,filesep,position,'_CFP_',num2str(frameNum),'.tif']));
    imRFPraw=double(imread([folder,filesep,position,'_Texas Red_',num2str(frameNum),'.tif']));
    
    jitterx=100; jittery=100;
    imCFPcrop=CropJitter(imCFPraw,jitterx,jitterx,jittery,jittery,0,0);
    imRFPcrop=CropJitter(imRFPraw,jitterx,jitterx,jittery,jittery,0,0);
       
    %%%%%% Background-subtract CFP/FRET images
    bgmaskCFP=getBGMask(imCFPcrop);
    imCFPbg=subBG(imCFPcrop,bgmaskCFP);
    bgmaskRFP=getBGMask(imRFPcrop);
    imRFPbg=subBG(imRFPcrop,bgmaskRFP);
    %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMask((imRFPcrop),2000);
    maskFinal{frameNum}=mask;
    cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Detrmine ratio
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',1),'replicate');
    imRFPbg(~mask)=nan;
    imRFP=ndnanfilter(imRFPbg,fspecial('disk',1),'replicate');
    imRatioTemp=imCFP./imRFP;
    imRatioTemp(~mask)=nan;
    imRatio_raw{frameNum}=imRatioTemp;
    %%%%%% Determine scaling for representation
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),3),1),round(prctile(imRatioTemp(:),97),1)];
    end
    %%%%%% Generate and write files for raw ratio and outlined objects
    tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imRFPOutline{frameNum}=DrawMaskOutline(imRFPcrop,mask);
    %imRFPOutlineTemp=imRFPOutline{frameNum};
    %imRFPOutline{frameNum}=DrawMaskOutlineRuby(imRFP_raw,mask);
    %imwrite(imRFPOutline,[datadir,filesep,position,'_FRET.tif'],'WriteMode','append','Compression','none');
    %imwrite(tempRATIO,[datadir,filesep,position,'_RATIO_raw_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
    %imwrite(imRFPOutlineTemp,[datadir,filesep,position,'_RFP.tif'],'WriteMode','append','Compression','none');
end

for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
end
save([datadir,filesep,position,'_RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','imRFPOutline');
save([datadir,filesep,position,'_Bleach_raw.mat'],'bleach_raw');

end
