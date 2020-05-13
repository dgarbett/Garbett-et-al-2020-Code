%% ParallelFRET_NG.m
% This script sets up parallel processing to obtain FRET data of time-lapse
% series by calling getFRETData in a parfor loop.
% For all movies in the same folder, a 
% Arnold, 9 June 2015

%% Parameter setup
clear;clc;close all; 

startpos=0;
endpos=6;
j=0;
for i=startpos:endpos
    j=j+1;
    position{j}=['Pos',num2str(i)];
end
%%
parfor k=1:length(position)
    
    getFRETDataNG_RFP(position{k});

end
disp('done!');

%% Bleaching correction
% calculate average bleaching from all movies in the same 
datadir='K:\data\150603_NG\RhoA2G_CDH5-mRuby3_timelapse_1min\data_150611_2';
close;
for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
%     bleach_raw(bleach_raw>3)=NaN;
%     bleach_raw(bleach_raw<1)=NaN;
    bleach_raw_all(:,k)=bleach_raw;
end
bleach_mean=nanmedian(bleach_raw_all,2);
plot(1:60,bleach_raw_all);hold on
plot(1:60,bleach_mean);

timepts=1:length(bleach_mean);
[fitvals]=polyfit(timepts',bleach_mean,1);
corr=polyval(fitvals,timepts);

plot(timepts,polyval(fitvals,timepts),'LineWidth',3);axis([0 length(bleach_mean) 1 2]);

%% correct bleaching
slope=fitvals(1);
parfor k=1:length(position)
    
   correctBleaching(position{k},slope);

end
disp('done!');



% %%%%%% Correct imRatio_raw for bleaching
% imRatio={};
% for frameNum=1:length(imRatio_raw)
%    if frameNum==1
%        frame_1=imRatio_raw{1};
%         colorRange = [round(prctile(frame_1(:),5),1),round(prctile(frame_1(:),97),1)];
% %        colorRange = [3.5, 5.5];
%    end
%    tempRATIO=imRatio_raw{frameNum}./corr(frameNum);
%    imRatio{frameNum}=tempRATIO;
%    tempRATIOforstack=ratio2RGB(tempRATIO,colorRange);%Cdc42
%    stitched=[imFRETOutline{frameNum} imRFPOutline{frameNum} tempRATIOforstack];
%    imwrite(stitched,[datadir,filesep,position,'_stitched_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
%    disp(num2str(frameNum));
% end       

