%% 

clear; clc;
root='K:\data\150701_IX';
rawdir=[root,filesep,'CR_thrombin'];
bgpath=[root,filesep,'backgroundCR'];
%imagepath=[root,filesep,'CDH5-Ruby_CC-Rac1_1min_4'];
imagepath=[root];
% %% Generate averaged background images

channels={'CFP' 'FRET'};
%channels={'CFP' 'DAPI' 'mCherry'};
%%
for chan=1:numel(channels);
    calculate_bg_img_rm_blobs(bgpath,channels{chan});
    disp(channels(chan));
    disp(num2str(chan));
end

%% Generate mean stack of first frame of all time-lapse series
CFPStack=[];
FRETStack=[]; k=0;
for row=2:3
    for col=8:9
        for site=1:4
            k=k+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            CFP_temp=imread([rawdir,filesep,shot,'_CFP_3.tif']);
            FRET_temp=imread([rawdir,filesep,shot,'_FRET_3.tif']);
            CFPStack(:,:,k)=CFP_temp;
            FRETStack(:,:,k)=FRET_temp;
            
        end
    end
end

CFP_AV=uint16(mean(CFPStack,3));
FRET_AV=uint16(mean(FRETStack,3));

imwrite(CFP_AV,[bgpath,filesep,'AVG_rawdata_CFP.tif'],'TIFF','Compression','None');
imwrite(FRET_AV,[bgpath,filesep,'AVG_rawdata_FRET.tif'],'TIFF','Compression','None');


%% Image alignment using Sean's method

alignStack(:,:,2)=imread([bgpath,filesep,'AVG_rawdata_FRET.tif']);
alignStack(:,:,1)=imread([bgpath,filesep,'AVG_rawdata_CFP.tif']);
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([bgpath,filesep,'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');