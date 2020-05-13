function correctBleachingExp2(position,fitpara,datadir)
% What it does
% Normalizes imRatio_raw for each frame by the median intensity per frame of the entire
% time series. 
% bleaching correction by division by linear bleaching function normalized
% by the median of its values. 

load([datadir,filesep,position,'_RatioData_raw.mat']);
load([datadir,filesep,position,'_Bleach_raw.mat']);
timepts=1:length(imRatio_raw);
corr=feval(fitpara,timepts);
%corr=polyval([slope 1],timepts);
corr_norm=corr./median(corr);
%plot(timepts,corr_norm);

normfact=nanmedian(bleach_raw);
%corr=polyval([slope 1],timepts);
%corr_norm=corr./median(corr);
%colorRange=[0.8 1.2];
%normalize RatioData_raw by median value
for frameNum=1:length(imRatio_raw)
    tempRATIO_raw=imRatio_raw{frameNum};
    tempRATIO_norm=tempRATIO_raw./normfact;
    tempRATIO_corr=tempRATIO_norm./corr_norm(frameNum);
     if frameNum==1
        colorRange = [round(prctile(tempRATIO_corr(:),1)*10)/10,round(prctile(tempRATIO_corr(:),98)*10)/10];
     end
%     colorRange=[0.8 1.1];
    imRatio{frameNum}=tempRATIO_corr;
    tempRATIOforstack=ratio2RGB(tempRATIO_corr,colorRange);%Cdc42
   %stitched=[imFRETOutline{frameNum} imRFPOutline{frameNum} tempRATIOforstack];
    stitched=[imRFPOutline{frameNum} tempRATIOforstack];
%   stitched=[imFRETOutline{frameNum} tempRATIOforstack];
   imwrite(stitched,[datadir,filesep,position,'_stitched_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
%    imwrite(imFRETOutline{frameNum},[datadir,filesep,position,'_FRETOutline.tif'],'WriteMode','append','Compression','none');
%    imwrite(imRFPOutline{frameNum},[datadir,filesep,position,'_RFPOutline.tif'],'WriteMode','append','Compression','none');
%    imwrite(tempRATIOforstack,[datadir,filesep,position,'_RATIO_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
   
   disp(num2str(frameNum));
end       
   save([datadir,filesep,position,'_RatioData.mat'],'-v7.3','maskFinal','cellCoors','imRatio');
end    
