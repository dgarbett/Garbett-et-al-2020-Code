% Image processing for/ratio imaging 
%
% Input
%   -TIFF sequence from 3I microscope which is named as "Capture #_... 
%   - 
%
% Output
%   - movie file with ratio images
%   - cell array with ratio images stored
%
% Background subtraction (subBG.m
%   - After generating a mask containing FRET-probe expressing cells, the
%   background level is determined from background regions and normalized 
%   by using the background image to account for uneven illumination 
%   (flatfield correction)
%
% Adapted from FRET code and getMyoDepData from Arnold Hayer
% Damien Garbett
%
%
%% USER INPUTS, check carefully before running
clear 
clc 
close all
root = 'G:\05.11.16\JustRatio';
specfolder = 'Cap1Crop2'
capnum='1'; %Capture number 
SF=0; %start frame number
EF=299; %end frame number
iRFPchannel=0; %set to Channel number for iRFP (in log file) %set to Ftractin
GFPchannel=1; %set to Channel number for GFP (in log file) %set to PLS3
%YFPchannel=2; %set to Channel number for YFP (in log file)

%% Working folder
rawdir=[root,filesep,specfolder];
bgdir='G:\05.11.16\JustRatio\Cap1Crop2\background';
datadir=[rawdir,filesep,'ImageData_',date];
mkdir(datadir);
jittersize=50;
rx1=jittersize;rx2=jittersize;ry1=jittersize;ry2=jittersize; 

%% Get subfolders and filenames
shot=['Capture ' capnum]'; %gets the capture number defined above in user inputs
filenames=getFilenames(rawdir); 
frame=1;
C1filesearchstring=['T',num2str(frame,'%3.3d'),'_C',num2str(iRFPchannel)]; %search for iRFP channel
C2filesearchstring=['T',num2str(frame,'%3.3d'),'_C',num2str(GFPchannel)]; %search for GFP channel
%C3filesearchstring=['T',num2str(frameNum,'%3.3d'),'_C',num2str(YFPchannel)]; %search for YFP channel
C1filename=char(filenames(boolRegExp(filenames,C1filesearchstring))); %gets filename for iRFP
C2filename=char(filenames(boolRegExp(filenames,C2filesearchstring))); %gets filename for GFP
%C3filename=char(filenames(boolRegExp(filenames,C3filesearchstring))); %gets filename for YFP
iRFP_raw=double(imread([rawdir,filesep,C1filename]));
GFP_raw=double(imread([rawdir,filesep,C2filename]));
%YFP_raw=double(imread([rawdir,filesep,C3filename]));
%% %%%% Call background images
%GFPbgraw=double(imread([bgdir,filesep,'AVG_bgGFP.tif']));
GFPbgraw=double(imread([bgdir,filesep,'AVG_bgPLS3.tif'])); %changed to PLS3
%YFPbgraw=double(imread([bgdir,filesep,'AVG_bgYFP.tif']));
%iRFPbgraw=double(imread([bgdir,filesep,'AVG_bgiRFP.tif']));
iRFPbgraw=double(imread([bgdir,filesep,'AVG_bgFTR.tif'])); %changed to FTractin
%% Test parameters using a single frame
close all
shot=['Capture ' capnum]';
%range=[0.2 3.0];
%folder = [rawdir,filesep,shot];
frameNum=1;

iRFP_raw=double(imread([rawdir,filesep,C1filename]));
GFP_raw=double(imread([rawdir,filesep,C2filename]));
%YFP_raw=double(imread([rawdir,filesep,C3filename]));

%Background-subtract images
bgmask=getBGMask(GFP_raw+iRFP_raw);
imGFPbg=subBG(GFP_raw,bgmask,GFPbgraw);
%imYFPbg=subBG(YFP_raw,bgmask,YFPbgraw);
iRFPbg=subBG(iRFP_raw,bgmask,iRFPbgraw);

%using Arnolds masking method for YFP channel (getCellMask.m) input image and min cell size
sum=iRFP_raw+GFP_raw;
mask=getCellMask((sum),4000);
edgemask=bwperim(mask);

%my mask method adapted from Min's tutorial
% logimiRFP=log(1+iRFP_raw);
% normlogimiRFP=mat2gray(logimiRFP);
% autothreshold1=graythresh(normlogimiRFP);
% automask1=normlogimiRFP>autothreshold1;
% areas=cell2mat(struct2cell(regionprops(automask1,'area')))';
% mask1clean=bwareaopen(automask1,1000);
% mask=imfill(mask1clean, 'holes'); % Fill holes in mask

%Detrmine ratio of GFP/iRFP (GCamp6s/iRFP) after background subtraction
GFP_raw2=imGFPbg; %sets to background subtracted CFP
GFP_raw2(~mask)=nan;
imGFP=ndnanfilter(GFP_raw2,fspecial('disk',2),'replicate'); %changed disk to 2 to blur less
iRFP_raw2=iRFPbg; %sets to background subtracted iRFP
iRFP_raw2(~mask)=nan;
imRFP=ndnanfilter(iRFP_raw2,fspecial('disk',2),'replicate'); %changed disk to 2 to blur less
imRatio1=imGFP./imRFP;
imRatio1(~mask)=nan;
colorRange1 = [round(prctile(imRatio1(:),3),1),round(prctile(imRatio1(:),99),1)];
tempRATIO1=ratio2parula(imRatio1,colorRange1);
% Plot
figure; 
subplot(2,2,1); imshow(GFP_raw,[]);
subplot(2,2,2); assessmask(iRFP_raw,mask);
subplot(2,2,3); imshow(mask); 
subplot(2,2,4); imshow(tempRATIO1);
%pause(20);

%% Make movie
close all
shot=[specfolder]';
%range=[0.2 3.0];
filenames=getFilenames(rawdir); 

imRatioStack=[];maskFinal={};cellCoors={};
tic;
for frameNum=SF:EF
    disp(num2str(frameNum));
    C1filesearchstring=['T',num2str(frameNum,'%3.3d'),'_C',num2str(iRFPchannel)]; %search for iRFP channel
    C2filesearchstring=['T',num2str(frameNum,'%3.3d'),'_C',num2str(GFPchannel)]; %search for GFP channel
    %C3filesearchstring=['T',num2str(frameNum,'%3.3d'),'_C',num2str(CFPchannel)]; %search for YFP channel
    C1filename=char(filenames(boolRegExp(filenames,C1filesearchstring))); %gets filename for mCherry
    C2filename=char(filenames(boolRegExp(filenames,C2filesearchstring))); %gets filename for CFP
    %C3filename=char(filenames(boolRegExp(filenames,C3filesearchstring))); %gets filename for YFP
    iRFP_raw=double(imread([rawdir,filesep,C1filename]));
    GFP_raw=double(imread([rawdir,filesep,C2filename]));
    %YFP_raw=double(imread([rawdir,filesep,C3filename]));
    disp(C1filename);
    disp(C2filename);
    %disp(C3filename);
    AdjFrameCount=frameNum+1; %for counting because 3i starts at frame 0
    
    %Background-subtract images
    bgmask=getBGMask(GFP_raw+iRFP_raw);
    imGFPbg=subBG(GFP_raw,bgmask,GFPbgraw);
    %imYFPbg=subBG(YFP_raw,bgmask,YFPbgraw);
    iRFPbg=subBG(iRFP_raw,bgmask,iRFPbgraw);
    
    %Arnold mask method
    sum=iRFP_raw+GFP_raw;
    mask=getCellMask((sum),4000);
    edgemask=bwperim(mask);
    
    
    %my mask method adapted from Min's tutorial
%     logimiRFP=log(1+iRFP_raw);
%     normlogimiRFP=mat2gray(logimiRFP);
%     autothreshold1=graythresh(normlogimiRFP);
%     automask1=normlogimiRFP>autothreshold1;
%     areas=cell2mat(struct2cell(regionprops(automask1,'area')))';
%     mask1clean=bwareaopen(automask1,1000);
%     mask=imfill(mask1clean, 'holes'); % Fill holes in mask
    
    %[mask cellCoorsTemp]=getCellMask((mCherry_raw),1500);
    %maskFinal{AdjFrameCount}=mask;
    %cellCoors{AdjFrameCount}=cellCoorsTemp;

    % Detrmine ratio 1
    GFP_raw2=imGFPbg; %sets to background subtracted CFP
    GFP_raw2(~mask)=nan;
    imGFP=ndnanfilter(GFP_raw2,fspecial('disk',2),'replicate'); %changed disk to 2 to blur less
    iRFP_raw2=iRFPbg; %sets to background subtracted iRFP
    iRFP_raw2(~mask)=nan;
    imRFP=ndnanfilter(iRFP_raw2,fspecial('disk',2),'replicate'); %changed disk to 2 to blur less
    imRatio=imGFP./imRFP;
    imRatio(~mask)=nan;

    imRatioScaled(~mask)=0;
    imRatioStack(:,:,AdjFrameCount)=imRatio;
    
    if frameNum==0
        colorRange = [round(prctile(imRatio(:),3),1),round(prctile(imRatio(:),99),1)];
    end

    tempRATIO=ratio2parula(imRatio,colorRange);
      
    imwrite(tempRATIO,[datadir,filesep,shot','_Ratio.tif'],'WriteMode','append','Compression','none');
    imwrite(edgemask,[datadir,filesep,shot','_mask.tif'],'WriteMode','append','Compression','none');
    %imwrite(%imwrite(tempCFP,[root filesep 'data' filesep position,'_DAPI.tif'],'WriteMode','append','Compression','none');
    %imwrite(tempRATIO,[datadir,filesep,shot,'_RATIO.tif'],'WriteMode','append','Compression','none');
end     
%for Num=1:length(imRatioStack)
%    bleach_raw(Num)=nanmean(vect(imRatioStack{Num}));
%end
    
save([datadir,filesep,shot','_ratio_raw.mat'],'maskFinal','cellCoors','imRatioStack','mask');
%save([datadir,filesep,shot,'_Bleach_raw.mat'],'bleach_raw');
toc;
disp('done!');

