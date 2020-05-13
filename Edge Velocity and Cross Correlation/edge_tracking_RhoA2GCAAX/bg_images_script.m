%% 

clear; clc;

%% Generate background images

bgpath='K:\data\150603_NG\background';

channels={'EPI-FRET' 'EPI-RFP'};
%channels={'CFP' 'DAPI' 'mCherry'};

for chan=1:numel(channels);
    calculate_bg_NKG(bgpath,channels{chan});
    disp(channels(chan));
    disp(num2str(chan));
end

%% Get alignment parameters. 
%% Determine alignment parameters for the two channels
filenames=getSubdirectories(alignFolder);
filenames=filenames(boolRegExp(filenames,'^img*'));
for i=1:length(filenames)
    alignIm(:,:,i)=imread([alignFolder filesep filenames{i}]);
end
alignIm=sum(alignIm,3);
alignIm=imread([alignFolder,filesep,'AVG_stack.tif']);
alignStack=dualview2stack(alignIm);
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
save([root filesep 'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');

%% Plot alignment parameters
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar
