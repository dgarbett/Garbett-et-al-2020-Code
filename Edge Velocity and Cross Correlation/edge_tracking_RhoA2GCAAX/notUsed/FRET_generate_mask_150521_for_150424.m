% Image processing for/ratio imaging 
%
% Input
%   - Movie sequence with one image per timepoint, acquired using
%   beamsplitter
%   - A background image acquired the same way, ideally an average of
%   multiple background images. 
%   - Images of the alignment grid to calculate the offset parameters
%
% Output
%   - movie file with ratio images
%   - cell array with ratio images stored
%
% Background subtraction (subBG.m)
%   - After generating a mask containing FRET-probe expressing cells, the
%   background level is determined from background regions and normalized 
%   by using the background image to account for uneven illumination 
%   (flatfield correction)
%
% Arnold Hayer, 150115
%
% Adapted functions for use of the Hamamatsu sCMOS camera. 


%% Working folder
clear 
clc 
close all
%warning off
root = 'K:\data\150424_IX\RhoFRET-CDH5-CHY-2min';
rawdir=[root,filesep,'rawdata'];
bgdir='K:\data\150424_IX\RhoFRET-CDH5-CHY-2min\background-bin1';
datadir=[root,filesep,'data_150521'];
mkdir(datadir);

%% Generate averaged image of first frame for alignment
CFP_stack=[];FRET_stack=[];i=0;
for row=6:7
    for col=8:10
        for site=1:2
            i=i+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            CFP_temp=imread([rawdir,filesep,shot,'_CFP_3.tif']);
            FRET_temp=imread([rawdir,filesep,shot,'_FRET_3.tif']);
            CFP_stack(:,:,i)=CFP_temp;
            FRET_stack(:,:,i)=FRET_temp;
        end 
    end
end
CFP_AV=uint16(mean(CFP_stack,3));
FRET_AV=uint16(mean(FRET_stack,3));
imwrite(CFP_AV,[bgdir,filesep,'AVG_rawdata_CFP.tif'],'TIFF','Compression','None');
imwrite(FRET_AV,[bgdir,filesep,'AVG_rawdata_FRET.tif'],'TIFF','Compression','None');

%% Image alignment using Sean's method
imCFP_raw=imread([bgdir,filesep,'AVG_rawdata_CFP.tif']);
imFRET_raw=imread([bgdir,filesep,'AVG_rawdata_FRET.tif']);
alignStack(:,:,1)=imCFP_raw; alignStack(:,:,2)=imFRET_raw;
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');

%% If previously determined, load alignement parameters, test parameters
load([root,filesep,'alignment parameters pX pY.mat'],'pX','pY');
frameNum=1;
%Get subfolders and filenames
shot='4_2_1';
CFP_files=dir([rawdir,filesep,shot,'_CFP_*']);
FRET_files=dir([rawdir,filesep,shot,'_FRET_*']);
imCFP_raw=double(imread([rawdir,filesep,shot,'_CFP_',num2str(frameNum),'.tif']));
imFRET_raw=double(imread([rawdir,filesep,shot,'_FRET_',num2str(frameNum),'.tif']));
imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
subplot(1,2,1), showImagesMergeChannels(imCFP_raw,imFRET_raw);
subplot(1,2,2), showImagesMergeChannels(imaligned(:,:,1),imaligned(:,:,2));

%% Load and process background image (an average of multiple background images)
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
CFPbg_raw=double(imread([bgdir,filesep,'AVG_bgCFP.tif']));
FRETbg_raw=double(imread([bgdir,filesep,'AVG_bgFRET.tif']));
bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);

%% Test parameters using a single frame
close all
shot= '4_2_1';
%range=[0.1 0.25];
folder = [rawdir,filesep,shot];
frameNum=1;
imCFP_raw=double(imread([rawdir,filesep,shot,'_CFP_',num2str(frameNum),'.tif']));
imFRET_raw=double(imread([rawdir,filesep,shot,'_FRET_',num2str(frameNum),'.tif']));
imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);

imCFP_raw=imaligned(:,:,1);
imFRET_raw=imaligned(:,:,2);

bgmask=getBGMask(imCFP_raw+imFRET_raw);
imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);

subplot(1,2,2);imagesc(imFRETbg+imCFPbg);
subplot(1,2,1);imagesc(bgmask);

imSum=(double(imCFP_raw+imFRET_raw));

% Get mask
mask=getCellMask((imFRET_raw),2000);
%mask_raw=segmentImageUsingThreshAndSeparate(imSum,'separate','','mincellarea',2000);
%mask=imfilter(mask_raw,fspecial('disk',2),'symmetric');

BWoutline = bwperim(mask);
Segout = mat2gray(imFRET_raw,[round(prctile(imFRET_raw(:),2),1),round(prctile(imFRET_raw(:),98),1)]);
Segout(BWoutline) = max(Segout(:));
figure, imshow(Segout);

% Detrmine ratio
imFRETbg(~mask)=nan;
imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
imCFPbg(~mask)=nan;
imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
imRatio=imFRET./imCFP;
imRatio(~mask)=nan;

imRatioScaled(~mask)=0;

colorRange = [round(prctile(imRatio(:),2),1),round(prctile(imRatio(:),98),1)];

tempRATIO=ratio2RGB(imRatio,colorRange);

% Plot
imshowpair(Segout,tempRATIO,'montage');

%% Make movie
load([root,filesep,'alignment parameters pX pY.mat'],'pX','pY');
for row=4%6:7
    for col=2%5:7
        for site=1%:2
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            disp([num2str(row),'_',num2str(col),'_',num2str(site)]);
            
            binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
            CFPbg_raw=double(imread([bgdir,filesep,'AVG_bgCFP.tif']));
            FRETbg_raw=double(imread([bgdir,filesep,'AVG_bgFRET.tif']));
            bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
            bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
            CFPbg=bg2(:,:,1);
            FRETbg=bg2(:,:,2);

            CFP_files=dir([rawdir,filesep,shot,'_CFP_*']);

            imRatioStack=[];
            for frameNum=1:5%length(CFP_files)
                disp(num2str(frameNum));
                imCFP_raw=double(imread([rawdir,filesep,shot,'_CFP_',num2str(frameNum),'.tif']));
                imFRET_raw=double(imread([rawdir,filesep,shot,'_FRET_',num2str(frameNum),'.tif']));

                imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
                imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);

                imCFP_raw=imaligned(:,:,1);
                imFRET_raw=imaligned(:,:,2);

                bgmask=getBGMask(imCFP_raw+imFRET_raw);
                imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
                imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);

                imCFPbg=subBG(imCFP_raw,bgmask);
                imFRETbg=subBG(imFRET_raw,bgmask);
               
                [mask cellCoorsTemp]=getCellMask(imFRET_raw,2000);
                maskFinal{frameNum}=mask;
                cellCoors{frameNum}=cellCoorsTemp;
                
                % Detrmine ratio
                imFRETbg(~mask)=nan;
                imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
                imCFPbg(~mask)=nan;
                imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
                imRatioTemp=imFRET./imCFP;
                imRatioTemp(~mask)=nan;
                imRatio_raw{frameNum}=imRatioTemp;
                if frameNum==1
                   colorRange = [round(prctile(imRatioTemp(:),2),1),round(prctile(imRatioTemp(:),98),1)];
                end
                imRatioScaled(~mask)=0;
                tempRATIO=ratio2RGB(imRatioTemp,colorRange);
                
                BWoutline = bwperim(mask);
                Segout = mat2gray(imFRET_raw,[round(prctile(imFRET_raw(:),2),1),round(prctile(imFRET_raw(:),98),1)]);
                Segout(BWoutline) = max(Segout(:));
                %figure, imshow(Segout);

                imwrite(Segout,[datadir,filesep,shot,'_FRET.tif'],'WriteMode','append','Compression','none');
                %imwrite(tempDAPI,[root filesep 'data' filesep position,'_DAPI.tif'],'WriteMode','append','Compression','none');
                %imwrite(tempRATIO,[datadir,filesep,shot,'_RATIO.tif'],'WriteMode','append','Compression','none');
                imwrite(tempRATIO,[datadir,filesep,shot,'_RATIO_raw_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
            end
            
            % Detrmine linear fit parameters for FRET/CFP decay
            bleach_1=nanmean(vect(imRatio_raw{1}));
            for frameNum=1:length(imRatio_raw)
               bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
            end
            bleach_corr=bleach_raw/bleach_1;
            timepts=1:length(imRatio_raw);
            [fitvals]=polyfit(timepts,bleach_corr,1);

            corr=polyval(fitvals,timepts);

            %plot(timepts,polyval(fitvals,timepts));
            imRatio={};
            for frameNum=1:length(imRatio_raw)
               if frameNum==1
                   frame_1=imRatio_raw{1};
                   colorRange = [round(prctile(frame_1(:),5),1),round(prctile(frame_1(:),95),1)];
               end
               imRatio{frameNum}=imRatio_raw{frameNum}./corr(frameNum);
               tempRATIO=ratio2RGB(imRatio{frameNum},colorRange);%Cdc42
               imwrite(tempRATIO,[datadir,filesep,shot,'_RATIO_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
               disp(num2str(frameNum));
            end       
                        
            save([datadir,filesep,'_RatioData_raw.mat'],'maskFinal', 'cellCoors','imRatio_raw', 'imRatio');
            
        end %site
    end %col
end %row
