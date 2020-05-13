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
root = '..';
alignFolder=[root filesep 'alignment40x'];
mkdir([root,filesep,'data']);
datadir=([root,filesep,'data']);
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

%% If previously determined, load alignement parameters
load([root filesep 'alignment parameters pX pY.mat'],'pX','pY');

%% Test alignment parameters on a test image

imraw=imread(['..\sequence\img_000000000_EPI-FRET (dual view)_000.tif']);
imstack=dualview2stack(imraw);
imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
subplot(1,2,1), showImagesMergeChannels(imstack(:,:,1),imstack(:,:,2));
subplot(1,2,2), showImagesMergeChannels(imaligned(:,:,1),imaligned(:,:,2));

%% Load and process background image (an average of multiple background images)
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
bgraw=double(imread([root filesep 'background40x' filesep 'AVG_Stack.tif']));
bgraw_smooth=imfilter(bgraw,fspecial('disk',10),'symmetric');
bg1=dualview2stack(bgraw_smooth);
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);
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
imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);

subplot(1,2,2);imagesc(imFRETbg+imCFPbg);
subplot(1,2,1);imagesc(bgmask);

imSum=(double(imCFP_raw+imFRET_raw));

mask_raw=segmentImageUsingThreshAndSeparate(imSum,'separate','','mincellarea',2000);
mask=imfilter(mask_raw,fspecial('disk',2),'symmetric');

% Detrmine ratio
imFRETbg(~mask)=nan;
imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
imCFPbg(~mask)=nan;
imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
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
    imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
    imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);

    subplot(1,2,2);imagesc(imFRETbg+imCFPbg);
    subplot(1,2,1);imagesc(bgmask);

    imSum=double(imCFP_raw+imFRET_raw);
    % Get mask
    mask_raw=segmentImageUsingThreshAndSeparate(imSum,'separate','','mincellarea',2000);
    mask=imfilter(mask_raw,fspecial('disk',2),'symmetric');
    % Detrmine ratio
    imFRETbg(~mask)=nan;
    imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
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