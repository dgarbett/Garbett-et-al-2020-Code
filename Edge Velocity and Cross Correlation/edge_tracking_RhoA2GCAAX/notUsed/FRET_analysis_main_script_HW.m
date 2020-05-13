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
% Arnold Hayer, 131122
%
% 140602 - adapted for NKG data taken using the EMCCD camera - CFP and FRET
% images are switched relative to the previous configuration.
%% Preparation
clear; clc; close all;
root = 'D:\Stanford\Matlab\XCorr\From Arnold For cross-correlation between FRET and edge dynamics (150109)\FRET_analysis';
addpath([root,filesep,'scripts/trackingcode']);
alignFolder=[root filesep 'alignment40x'];
mkdir([root filesep,'data']);
datadir=([root filesep,'data']);
binning=1; % needs to be changed only when alignment images and data were acquired at distinct binning seetings

jittersize=10;
rx1=jittersize;rx2=jittersize;ry1=jittersize;ry2=jittersize;
% %% Determine alignment parameters for the two channels
% % filenames=getFilenames(alignFolder);
% % filenames=filenames(boolRegExp(filenames,'^img*'));
% % for i=1:length(filenames)
% %     alignIm(:,:,i)=imread([alignFolder filesep filenames{i}]);
% % end
% % alignIm=sum(alignIm,3);
% alignIm=imread('M:\data\141121_NG\test\AVG_Stack.tif');
% alignStack=dualview2stack(alignIm);
% [pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
% save([root filesep 'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');

% %% Plot alignment parameters
% figure;
% subplot(1,2,1); imagesc(dxMat1); colorbar
% subplot(1,2,2); imagesc(dyMat1); colorbar
%% Parameters for cell edge parametrization
nFretWindows=100;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=3;                                      % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge
edgeDepthDist=20;                                   % Number of pixels deep for the windows for computing FRET values.
dataFolder=[root filesep 'Data'];
%% If previously determined, load alignement parameters
load([root filesep 'alignment parameters pX pY.mat'],'pX','pY');
%% Test alignment parameters on a test image
imraw=imread([root '\sequence\img_000000000_EPI-FRET (dual view)_000.tif']);
imstack=dualview2stack(imraw);
imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
subplot(1,2,1), showImagesMergeChannels(imstack(:,:,1),imstack(:,:,2));
title('before alignment');
subplot(1,2,2), showImagesMergeChannels(imaligned(:,:,1),imaligned(:,:,2));
title('after alignment');
%% Load and process background image (an average of multiple background images)
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
bgraw=double(imread([root filesep 'background40x' filesep 'AVG_Stack.tif']));
bgraw_smooth=imfilter(bgraw,fspecial('disk',10),'symmetric');
bg1=dualview2stack(bgraw_smooth);
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);
%% Compute masks and background and cell coordinates
subfolder=[root filesep 'sequence'];
fprintf('%s\n',subfolder);
filenames=getFilenames(subfolder,'.tif$');
% if ~exist([subfolder filesep 'masks and coors.mat'],'file')
%load([subfolder filesep 'adjusted alignment parameters pX pY.mat'],'pX','pY');
clear maskFinal; clear coors; clear maskBGinitial
vw=VideoWriter([subfolder filesep 'movie.avi']);
vw.set('FrameRate',8);
vw.set('Quality',100);
open(vw); figure('Position',[50 50 600 700]);
for imnum=1:length(filenames)
    im1=imread([subfolder filesep filenames{imnum}]);
    im1s=dualview2stack(im1);
    im1s=dualviewAlignFromFittedSurface(im1s,pX,pY,binning);
    imCFP_raw=im1s(:,:,1);
    imFRET_raw=im1s(:,:,2);
    imSum=sum(double(im1s),3);
    bgmask{imnum}=getBGMask(imSum);
    imCFP=subBG(imCFP_raw,bgmask{imnum},CFPbg);
    imFRET=subBG(imFRET_raw,bgmask{imnum},FRETbg);
    imCFP=ndnanfilter(imCFP,fspecial('gaussian',7,3),7);
    imFRET=ndnanfilter(imFRET,fspecial('gaussian',7,3/binning),7);
    %         binning=round(1040/size(imSum,1));
    mask_raw=segmentImageUsingThreshAndSeparate(imSum,'separate','','mincellarea',2000);
    [maskFinal{imnum},coors{imnum}]=fretGetCellMasks(2000,imCFP,imFRET,mask_raw); %% gets FRET and mask values 2000 min cell size
    imCFP(~maskFinal{imnum})=nan;
    imFRET(~maskFinal{imnum})=nan;
    %         if binning==1
    %             imFRET=imageBin(imFRET,2); imCFP=imageBin(imCFP,2);
    %         end
    imRatio=imFRET./imCFP;
    if imnum==1
        figure('Position',[100 100 size(imRatio,2) size(imRatio,1)-200]*1.1);
        colorRange = [1.4 2.4];%round(prctile(imRatio(:),1),1),round(prctile(imRatio(:),95),1)];
    end
    imRatio(imRatio<(colorRange(1)-0.1))=(colorRange(1)-0.1); %Distinguish background from low FRET cell
    imagesc(imRatio,[(colorRange(1)-0.1) (colorRange(2)+0.1)]); colormap('jet_black'); colorbar;
    title(['Cdc42 activity']);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    axis image
    writeVideo(vw,getframe(gcf));
end
close(vw);
save([subfolder filesep 'masks and coors.mat'],'bgmask','maskFinal','coors');
% end

%% Compute and Save Cell Edge Dynamics Data
load([root filesep 'sequence' filesep 'masks and coors.mat']);
numBefore=0; % e.g. frames before stimulation
clear edgeData;
cellCount=0;
timeinterval=15;
%     expFolders=getSubdirectories(root);
%         for snum=1:length(expFolders)
%             subfolder=[folder filesep subdir{snum}];
subfolder=[root filesep 'sequence'];
fprintf('%s\n',subfolder);
filenames=getFilenames(subfolder,'.tif$'); filenames=filenames(~boolRegExp(filenames,'RFP'));
try
    timeData=fretReadTimeDataFromMetadata(subfolder,numBefore);
catch
    timeData=(0:timeinterval:timeinterval*(length(filenames)-1))'; %% AH> 
end
%load([subfolder filesep 'masks and coors.mat'],'maskBGinitial','maskFinal','coors');

load([subfolder filesep 'masks and coors.mat'],'maskFinal','coors');

% perform tracking of mask centroid
if min(cellfun(@(x) size(x,1),coors))>0
    traj=ultTrackAnnSearch(coors,'pairrule','fwdbckmtch','maxdisp',100,'verbose',false);
else
    traj={};
end

fprintf('%i trajectories.\n',length(traj));
for trajNum=1:length(traj)
    if size(traj{trajNum},1)>=20
        cellCount=cellCount+1;
        frameCount=0;
        clear edgeCoors; clear edgeCoorsSmoothed; clear windowCoors; clear protvals; clear fretvals; clear indFrontCenter; clear frontMap; clear backMap; clear windowCoors
        startFrame=traj{trajNum}(1,end);
        endFrame=traj{trajNum}(end,end);
        for imnum=startFrame:endFrame
            frameCount=frameCount+1;
            
            % Read images and compute FRET values
            im1=imread([subfolder filesep filenames{imnum}]);
            im1s=dualview2stack(im1);
            im1s=dualviewAlignFromFittedSurface(im1s,pX,pY,binning);
            %                         binning=round(1040/size(im1s,1));
            imCFP=subBG(im1s(:,:,2),bgmask{imnum},CFPbg);
            imFRET=subBG(im1s(:,:,1),bgmask{imnum},FRETbg);
            
            % Identify the correct cell and make a one cell mask
            objects=regionprops(maskFinal{imnum},'PixelIdxList','PixelList','Centroid','BoundingBox');
            thisTraj=traj{trajNum};
            cellCent=round(thisTraj(find(thisTraj(:,end)==imnum,1),1:2));
            cellNum=find(round(arrayfun(@(x) x.Centroid(1),objects))==cellCent(1) & round(arrayfun(@(x) x.Centroid(2),objects))==cellCent(2));
            thisMask=false(size(maskFinal{imnum}));
            thisMask(objects(cellNum).PixelIdxList)=true;
            
            % Parametrize cell edge and compute protrusion values
            if frameCount==1
                [edgeCoors{frameCount},edgeCoorsSmoothed{frameCount}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning));
            else
                [edgeCoors{frameCount},edgeCoorsSmoothed{frameCount}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning),edgeCoors{frameCount-1});
                windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                protvals(:,frameCount-1)=vect(computeProtrusionValues(edgeCoors{frameCount-1},edgeCoorsSmoothed{frameCount-1},edgeCoors{frameCount}));
            end
            windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
            
            % Compute edge regions and local fret ratio values
            labelMask=getWindowLabelMap(thisMask,windowCoors{frameCount},edgeDepthDist/binning);
            for k=1:size(windowCoors{frameCount},1)
                fretvals(k,frameCount)=sum(imFRET(labelMask==k))/sum(imCFP(labelMask==k));
            end
        end
        protvalsWindow=zeros(size(fretvals)-[0 1]);
        for k=1:edgeOversamplingParam
            protvalsWindow=protvalsWindow+protvals(k:edgeOversamplingParam:end,:);
        end
        % Normalize protrusion data by dividing by the edgeOversampling factor and dividing by the time between images.
        protvalsWindow=(protvalsWindow/edgeOversamplingParam)./repmat(diff(timeData(startFrame:endFrame))',[nFretWindows 1]);
        
        % Fit a bleaching parameter
        bleachCorrectionLinear=fitline(1:size(fretvals,2),nanmean(fretvals,1));
        bleachCorrectionExp=fitline(1:size(fretvals,2),log(nanmean(fretvals,1)));
        
        [indFrontCenter,frontMap,backMap]=computeCenterOfCellFront(protvalsWindow);
        
        % Build an edgeData structure
        edgeData(cellCount).edgeCoors=edgeCoors;
        edgeData(cellCount).edgeCoorsSmoothed=edgeCoorsSmoothed;
        edgeData(cellCount).windowCoors=windowCoors;
        edgeData(cellCount).cellArea=traj{trajNum}(:,3)*binning*binning;
        edgeData(cellCount).protvalsAll=protvals;
        edgeData(cellCount).protvalsWindow=protvalsWindow;
        edgeData(cellCount).fretvals=fretvals;
        edgeData(cellCount).indFrontCenter=indFrontCenter;
        edgeData(cellCount).frontMap=frontMap;
        edgeData(cellCount).backMap=backMap;
        edgeData(cellCount).rootDir=root;
        edgeData(cellCount).folder=subfolder;
        edgeData(cellCount).trajNum=trajNum;
        edgeData(cellCount).startFrame=startFrame;
        edgeData(cellCount).endFrame=endFrame;
        edgeData(cellCount).bleachCorrectionLinear=bleachCorrectionLinear;
        edgeData(cellCount).bleachCorrectionExp=bleachCorrectionExp;
        bleachCorrAllValues.m(cellCount,1)=edgeData(cellCount).bleachCorrectionLinear.m;
        bleachCorrAllValues.b(cellCount,1)=edgeData(cellCount).bleachCorrectionLinear.b;
    end
end
% Compare bleaching correction values
bleachCorrAllValues.rel=bleachCorrAllValues.m./bleachCorrAllValues.b;
% save([dataFolder filesep '.mat'],'edgeData','bleachCorrAllValues');  % Save the data for this experimental group
%     end


















%% Choose subfolder, filnames, range
filenamecontains='^*img*';

% Test parameters using a single frame
close all
frameNum=1;
range=[1 3];
folder=[root filesep 'sequence'];
filenames=getFilenames(folder);
filenames=filenames(boolRegExp(filenames,filenamecontains));
disp(num2str(frameNum));
% Load image
imraw=double(imread([folder filesep filenames{frameNum}]));
im1=dualview2stack(imraw);
showImagesMergeChannels(im1(:,:,1),im1(:,:,2));
im2=dualviewAlignFromFittedSurface(im1,pX,pY,binning);
showImagesMergeChannels(im2(:,:,1),im2(:,:,2));
imCFP_raw=im2(:,:,1);
imFRET_raw=im2(:,:,2);

bgmask=getBGMask(imCFP_raw+imFRET_raw);
imCFP=subBG(imCFP_raw,bgmask,CFPbg);
imFRET=subBG(imFRET_raw,bgmask,FRETbg);

subplot(1,2,2);imagesc(imFRET+imCFP);
subplot(1,2,1);imagesc(bgmask);

imSum=(double(imCFP_raw+imFRET_raw));

mask_raw=segmentImageUsingThreshAndSeparate(imSum,'separate','','mincellarea',2000);
mask=imfilter(mask_raw,fspecial('disk',2),'symmetric');

% Detrmine ratio
imFRET(~mask)=nan;
imFRET=ndnanfilter(imFRET,fspecial('disk',3),'replicate');
imCFP(~mask)=nan;
imCFP=ndnanfilter(imCFP,fspecial('disk',3),'replicate');
imRatio=imFRET./imCFP;
imRatio(~mask)=nan;

imRatioScaled(~mask)=0;
tempRATIO=ratio2RGB(imRatio,range);

% Plot
figure; subplot(2,2,1); imshow(mat2gray(imSum));
subplot(2,2,2); imshow(mask);
subplot(2,2,3); imshow(tempRATIO);
subplot(2,2,4); imshow(mat2gray(imSum.*mask));

%% make movie

% Choose subfolder, filnames, range
binning=1;
% allSubdir=getSubdirectories(root); %[root filesep 'CDH5-TS']);
% dataSubdir=allSubdir(boolRegExp(allSubdir, '^*CC-Cdc42*')) % displays the subfolders used for analysis

filenamecontains='^*img*';

% subNum=2;
position=0;
range=[1 3];
close all;
%background image
bgraw=double(imread([root filesep 'background40x' filesep 'AVG_Stack.tif']));
bg1=dualview2stack(bgraw);
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);

folder=[root,filesep,'sequence'];
filenames=getFilenames(folder);
filenames=filenames(boolRegExp(filenames,filenamecontains));

for frameNum=1:(length(filenames));
    disp(num2str(frameNum));
    imraw=double(imread([folder filesep filenames{frameNum}]));
    im1=dualview2stack(imraw);
    
    im2=dualviewAlignFromFittedSurface(im1,pX,pY,binning);
    
    imCFP_raw=im2(:,:,1);
    imFRET_raw=im2(:,:,2);
    
    bgmask=getBGMask(imCFP_raw+imFRET_raw);
    imCFP=subBG(imCFP_raw,bgmask,CFPbg);
    imFRET=subBG(imFRET_raw,bgmask,FRETbg);
    
    subplot(1,2,2);imagesc(imFRET+imCFP);
    subplot(1,2,1);imagesc(bgmask);
    
    imSum=double(imCFP_raw+imFRET_raw);
    % Get mask
    mask_raw=segmentImageUsingThreshAndSeparate(imSum,'separate','','mincellarea',2000);
    mask=imfilter(mask_raw,fspecial('disk',2),'symmetric');
    % Detrmine ratio
    imFRET(~mask)=nan;
    imFRET=ndnanfilter(imFRET,fspecial('disk',3),'replicate');
    imCFP(~mask)=nan;
    imCFP=ndnanfilter(imCFP,fspecial('disk',3),'replicate');
    imRatio=imFRET./imCFP;
    imRatio(~mask)=nan;
    
    imRatioScaled(~mask)=0;
    tempRATIO=ratio2RGB(imRatio,range);
    
    close;
    imRatioStack(:,:,frameNum)=imRatio;
    tempFRET=imadjust(mat2gray(imFRET_raw));
    
    tempRATIO=ratio2RGB(imRatio,range);%Cdc42
    
    imwrite(tempFRET,[datadir filesep '_FRET.tif'],'WriteMode','append','Compression','none');
    %imwrite(tempDAPI,[root filesep 'data' filesep position,'_DAPI.tif'],'WriteMode','append','Compression','none');
    imwrite(tempRATIO,[datadir filesep '_RATIO.tif'],'WriteMode','append','Compression','none');
    
    
end
disp('done!');