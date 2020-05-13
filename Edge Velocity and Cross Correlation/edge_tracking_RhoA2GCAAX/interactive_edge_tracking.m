function interactive_edge_tracking(position,bgdir,rawdir,datadir,jitter) %#ok<INUSD>
%%%%meant for fret images with the beam splitter
%run bg_images_alignment_script before hand... try to build in before hand
%% Set up
filesep= '/';
load([bgdir,filesep,'alignment parameters pX pY.mat']);
jitter=0;
jitterx=25;
jittery=60;
binning=1;
%position = 'Pos1';

%open images and perform the dual view alignment
folder=[rawdir,filesep,position];
filenames=getFilenames(folder);
filenamesFRET=filenames(boolRegExp(filenames,'FRET'));
filenamesRFP=filenames(boolRegExp(filenames,'RFP'));
imRatio_raw={};maskFinal={};cellCoors={};

load([bgdir,filesep,'alignment parameters pX pY.mat']);


%%for each position(k) open up the first image and last image
imrawRFP=double(imread([folder,filesep,filenamesRFP{1}]));
imRFP1=dualview2stack(imrawRFP);
imRFP2=dualviewAlignFromFittedSurface(imRFP1,pX,pY,binning);
imRFP_raw=imRFP2(:,:,1);

imrawRFP=double(imread([folder,filesep,filenamesRFP{end}]));
imRFP1=dualview2stack(imrawRFP);
imRFP2=dualviewAlignFromFittedSurface(imRFP1,pX,pY,binning);
imRFP_raw_end=imRFP2(:,:,1);
clear imRFP1 imRFP2;

figure(1)
h=imagesc(log(imRFP_raw));
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
for frameNum=1:length(filenamesFRET)
    
    %load in images
    imrawFRET=double(imread([folder,filesep,filenamesFRET{frameNum}]));
    imFRET1=dualview2stack(imrawFRET);
    imFRET2=dualviewAlignFromFittedSurface(imFRET1,pX,pY,binning);
    
    imrawRFP=double(imread([folder,filesep,filenamesRFP{frameNum}]));
    imRFP1=dualview2stack(imrawRFP);
    imRFP2=dualviewAlignFromFittedSurface(imRFP1,pX,pY,binning);
    
    imCFP_raw=imFRET2(:,:,2);
    imFRET_raw=imFRET2(:,:,1);
    imRFP_raw=imRFP2(:,:,1);
    
    % jitter correction using RFP image
    if jitter==1;
        if frameNum>1
            [alignpara ~]=dftregistration(fft2(imRFP_raw0),fft2(imRFP_raw),100);
            xx=alignpara(4);yy=alignpara(3);
            x=x+xx;y=y+yy;
        else
            x=0;y=0;    %store jitters data
        end
        imRFP_raw0=imRFP_raw;
        disp([x,y])
        
        imCFP_raw=CropJitter(imCFP_raw,jitterx,jitterx,jittery,jittery,x,y);
        imFRET_raw=CropJitter(imFRET_raw,jitterx,jitterx,jittery,jittery,x,y);
        imRFP_raw=CropJitter(imRFP_raw,jitterx,jitterx,jittery,jittery,x,y);
        
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
    
    maskAll = zeros(size(imRFP_raw));
    cellCoorsAll = [];
    %iterate through the cells to get data for each frame
    for c = 1:numCells;
        % First mask using RFP channel
        [mask cellCoorsTemp]=getIndividualCellMask((imRFP_raw.*cellMask{c}),2000);
        maskAll  = maskAll  | mask;
        cellCoorsAll = [cellCoorsAll; cellCoorsTemp];
    end %through cells
    imFRET_raw(~maskAll)=nan;
    imFRET=ndnanfilter(imFRET_raw,fspecial('disk',1),'replicate');
    imCFP_raw(~maskAll)=nan;
    imCFP=ndnanfilter(imCFP_raw,fspecial('disk',1),'replicate');
    imRatioTemp=imFRET./imCFP;
    
    
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),3)*10)/10,round(prctile(imRatioTemp(:),97)*10)/10];
    end
    tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,maskAll);
    imRFPOutline{frameNum}=DrawMaskOutlineRuby(imRFP_raw,maskAll);
    
    maskFinal{frameNum}=maskAll;
    cellCoors{frameNum}=cellCoorsAll;
    imRatio_raw{frameNum}=imRatioTemp;
end%through frames

for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
end

save([datadir,filesep,position,'_RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','imFRETOutline','imRFPOutline');
save([datadir,filesep,position,'_Bleach_raw.mat'],'bleach_raw');

end

