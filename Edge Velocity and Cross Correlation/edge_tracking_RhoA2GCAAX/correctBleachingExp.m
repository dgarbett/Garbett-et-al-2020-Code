function correctBleachingExp( position,slope,datadir)

% datadir=['K:\data\150617_NG\rawdata\RhoA2G_20s_2\data_150618'];
% load([datadir,filesep,position,'_RatioData_raw.mat']);
% timepts=1:length(imRatio_raw);
corr=polyval([slope 1],timepts);
% plot(timepts,corr);
%%%%%% Correct imRatio_raw for bleaching
imRatio={};
for frameNum=1:length(imRatio_raw)
   if frameNum==1
       frame_1=imRatio_raw{1};
        %colorRange = [round(prctile(frame_1(:),2),1),round(prctile(frame_1(:),98),1)];
        colorRange = [1.6 3];
   end
   tempRATIO=imRatio_raw{frameNum}./corr(frameNum);
   imRatio{frameNum}=tempRATIO;
   tempRATIOforstack=ratio2RGB(tempRATIO,colorRange);%Cdc42
   stitched=[imFRETOutline{frameNum} imRFPOutline{frameNum} tempRATIOforstack];
   %stitched=[imFRETOutline{frameNum} tempRATIOforstack];
   imwrite(stitched,[datadir,filesep,position,'_stitched_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
   disp(num2str(frameNum));
end       
   save([datadir,filesep,position,'_RatioData.mat'],'maskFinal','cellCoors','imRatio');

end

