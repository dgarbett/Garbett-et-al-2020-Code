function getFRETDataNG_RFP_junction(position,bgdir,rawdir,datadir,jitter)
% image processing for FRET data analysis 
%position='Pos1';
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
% root = 'K:\data\150617_NG';
% rawdir=[root,filesep,'rawdata',filesep,'RhoA2G_LynFRB-P2A-ARHGEF1-GEF_30s_1'];
% bgdir=[root,filesep,'background'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
% datadir=[root,filesep,'data_150617'];
% if ~exist(datadir)
%     mkdir(datadir);
% end
jitterx=25;
jittery=60;

%% %%%% Call background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
bgrawFRET=double(imread([bgdir,filesep,'AVG_bgEPI-FRET.tif']));
bgFRET1=dualview2stack(bgrawFRET);
bgFRET2=dualviewAlignFromFittedSurface(bgFRET1,pX,pY,binning);
CFPbg_raw=bgFRET2(:,:,2);
FRETbg_raw=bgFRET2(:,:,1);

bgrawRFP=double(imread([bgdir,filesep,'AVG_bgEPI-RFP.tif']));
bgRFP1=dualview2stack(bgrawRFP);
bgRFP2=dualviewAlignFromFittedSurface(bgRFP1,pX,pY,binning);
RFPbg_raw=bgRFP2(:,:,1);

%%%%%% Loop through frames
folder=[rawdir,filesep,position];
%folder=[rawdir,filesep,position,filesep,'Pos0'];
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
    %%%%%% Jitter correction (optional)
    if jitter==1;
        if frameNum>1
            [alignpara ~]=dftregistration(fft2(imRFP_raw0),fft2(imRFP_raw),100);  
            xx=alignpara(4);yy=alignpara(3);
%             if (abs(xx)+abs(yy)>jitterx)
%                 xx=0;yy=0;
%                 imRFP_raw=imRFP_raw0;
%             end
            x=x+xx;y=y+yy;
        else
            x=0;y=0;    %store jitters data
        end
        imRFP_raw0=imRFP_raw;
        disp([x,y])
        
        imCFP_raw=CropJitter(imCFP_raw,jitterx,jitterx,jittery,jittery,x,y);
        imFRET_raw=CropJitter(imFRET_raw,jitterx,jitterx,jittery,jittery,x,y);
        imRFP_raw=CropJitter(imRFP_raw,jitterx,jitterx,jittery,jittery,x,y);
        
        CFPbg=CropJitter(CFPbg_raw,jitterx,jitterx,jittery,jittery,x,y);
        FRETbg=CropJitter(FRETbg_raw,jitterx,jitterx,jittery,jittery,x,y);
        RFPbg=CropJitter(RFPbg_raw,jitterx,jitterx,jittery,jittery,x,y);
    else
        CFPbg=CFPbg_raw;
        FRETbg=FRETbg_raw;
        RFPbg=RFPbg_raw;
    end
        
    %%%%%% Background-subtract CFP/FRET images
    %bgmask=getBGMask(imCFP_raw);
    %imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
    %imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
%    keyboard
    %imRFPbg=subBG(imRFP_raw,bgmask,RFPbg);
    %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMask((imRFP_raw),2000);
    
    %erode image and get junction area of interest
    maskIn = imerode(mask,strel('disk',3,8));
    maskRing = mask;
    maskRing(maskIn) = 0;
    
    maskFinal{frameNum}=maskRing;
    cellCoors{frameNum}=cellCoorsTemp;
    
    
    
    %%%%%% Detrmine ratio
    imFRET_raw(~maskRing)=nan;
    imFRET=ndnanfilter(imFRET_raw,fspecial('disk',1),'replicate');
    imCFP_raw(~maskRing)=nan;
    imCFP=ndnanfilter(imCFP_raw,fspecial('disk',1),'replicate');
    imRatioTemp=imFRET./imCFP;
    %imRatioTemp(~mask)=nan;
    imRatio_raw{frameNum}=imRatioTemp;
    %%%%%% Determine scaling for representation
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),3)*10)/10,round(prctile(imRatioTemp(:),97)*10)/10];
    end
    %%%%%% Generate and write files for raw ratio and outlined objects
    tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);
    imRFPOutline{frameNum}=DrawMaskOutlineRuby(imRFP_raw,mask);
    %imwrite(imFRETOutline,[datadir,filesep,position,'_FRET.tif'],'WriteMode','append','Compression','none');
    %imwrite(tempRATIO,[datadir,filesep,position,'_RATIO_raw_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
    %imwrite(imRFPOutline,[datadir,filesep,position,'_RFP.tif'],'WriteMode','append','Compression','none');
end

for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
end
save([datadir,filesep,position,'_RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','imFRETOutline','imRFPOutline');
save([datadir,filesep,position,'_Bleach_raw.mat'],'bleach_raw');

% %% %%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay
% bleach_1=nanmean(vect(imRatio_raw{1}));
% for frameNum=1:length(imRatio_raw)
%    bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
% end
% bleach_corr=bleach_raw/bleach_1;
% timepts=1:length(imRatio_raw);
% [fitvals]=polyfit(timepts,bleach_corr,1);
% corr=polyval(fitvals,timepts);
% %plot(timepts,polyval(fitvals,timepts));
% %%%%%% Correct imRatio_raw for bleaching
% imRatio={};
% for frameNum=1:length(imRatio_raw)
%    if frameNum==1
%        frame_1=imRatio_raw{1};
%         colorRange = [round(prctile(frame_1(:),5),1),round(prctile(frame_1(:),97),1)];
% %        colorRange = [3.5, 5.5];
%    end
%    tempRATIO=imRatio_raw{frameNum}./corr(frameNum);
%    imRatio{frameNum}=tempRATIO;
%    tempRATIOforstack=ratio2RGB(tempRATIO,colorRange);%Cdc42
%    stitched=[imFRETOutline{frameNum} imRFPOutline{frameNum} tempRATIOforstack];
%    imwrite(stitched,[datadir,filesep,position,'_stitched_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
%    disp(num2str(frameNum));
% end       
% %%%%%% Save data
% save([datadir,filesep,position,'_RatioData.mat'],'maskFinal','cellCoors','imRatio');
 end

