%% For experiment 150424
% Wells 6_2, 7_2, 6_5, 7_5 are controls, 

%% Working folder
clear;clc;close all; 
root = 'K:\data\150701_IX';
rawdir=[root,filesep,'CY_thrombin'];
bgdir=[root,filesep,'backgroundCY'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[root,filesep,'data_150701'];

%% Plot average intensity per frame
i=0;j=0;
for row=2:3
    for col=8:9
        j=j+1;
        for site=1:3
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            i=i+1;
            load([datadir,filesep,shot,'_RatioData.mat']);
            %load([datadir,filesep,shot,'_Bleach_raw.mat']);
            %meanFRET(:,site)=bleach_raw;
            for frame=1:length(imRatio)
                %meanFRET(frame,site)=nanmedian(vect(imRatio{frame}));
                meanFRET(frame,site)=nanmedian(vect(imRatio{frame}));
            end
            disp(shot);
        end %site
        FRET_divide=mean(meanFRET(1:10,:),1);
        FRET_divide_mat=repmat(FRET_divide,70,1);
        meanFRET_norm=meanFRET./FRET_divide_mat;
        subplot(2,8,j);plot(1:70,meanFRET_norm); axis([0 70 0.9 1.3]);
    end %col
end %row
disp('done!');

% %% calculate mean FRET norm (normalized by mean of the first 10 timpoints)
% FRET_divide=mean(meanFRET(1:10,:),1);
% FRET_divide_mat=repmat(FRET_divide,60,1);
% 
% meanFRET_norm=meanFRET./FRET_divide_mat;
% load('timestamp_cols234.mat');
% timestamp234=timestmp;
% load('timestamp_cols567.mat');
% timestamp567=timestmp;
% i=0;
% for col=2:4
%    i=i+1;
%    subplot(2,6,i);
%    plot(timestamp234,mean
%    
%         
%     
% 
% 
% 
% 
% load('K:\data\150506_IXM\predrug_drug\data_150506\timestamp_RhoA.mat');
% timeaxis=timestmp(:,1);
% 
% subplot(1,3,1);plot(timeaxis,meanFRET_norm(:,1:2));axis([1 3000 0.75 1.3]);
% subplot(1,3,2);plot(timeaxis,meanFRET_norm(:,3:4));axis([1 3000 0.75 1.3]);
% subplot(1,3,3);plot(timeaxis,meanFRET_norm(:,5:6));axis([1 3000 0.75 1.3]);
% 
% 
