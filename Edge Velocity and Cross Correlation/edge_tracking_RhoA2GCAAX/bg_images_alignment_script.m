%% 

clear; clc;
root = '/Volumes/BisackUp/Meyer_Lab/CDH5-TS-Thrombin';
bgpath=[root,filesep,'background'];
imagepath=[root,filesep,'CDH5-TS-1min_1'];
%imagepath=[root];
%% Generate averaged background images

channels={'EPI-RFP' 'EPI-FRET' 'RFP'};
%channels={'CFP' 'DAPI' 'mCherry'};

for chan=1:numel(channels);
    calculate_bg_NKG(bgpath,channels{chan});
    disp(channels(chan));
    disp(num2str(chan));
end

%% Generate mean stack of first frame of all time-lapse series

FRETStack=[]; i=0;
subdir=getSubdirectories(imagepath);
subdir=subdir(boolRegExp(subdir,'Pos'));
for i=1:length(subdir)
    FRET_temp=imread([imagepath,filesep,subdir{i},filesep,'img_000000000_EPI-FRET (dual view)_000.tif']);
    FRETStack(:,:,i)=FRET_temp;
    
end

FRET_AV=uint16(mean(FRETStack,3));
imwrite(FRET_AV,[bgpath,filesep,'AVG_rawdata_FRET.tif'],'TIFF','Compression','None');


%% Image alignment using Sean's method
FRET_AV=imread([bgpath,filesep,'AVG_rawdata_FRET.tif']);
alignStack=dualview2stack(FRET_AV);
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([bgpath,filesep,'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');