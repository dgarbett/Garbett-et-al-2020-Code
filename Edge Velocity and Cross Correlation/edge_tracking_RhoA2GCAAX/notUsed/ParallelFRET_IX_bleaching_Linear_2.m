%% ParallelFRET_NG.m
% This script sets up parallel processing to obtain FRET data of time-lapse
% series by calling getFRETData in a parfor loop.
% For all movies in the same folder, a 
% Arnold, 9 June 2015

%% Parameter setup
clear;clc;close all; 
root = 'K:\data\150701_IX';
rawdir=[root,filesep,'CY_thrombin'];
bgdir=[root,filesep,'backgroundCY'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[root,filesep,'data_150702'];
if ~exist(datadir)
    mkdir(datadir);
end

%jitter=1; %jitter correction, 0 or 1.
k=0;
for row=2:3
    for col=2
        for site=1:3
            k=k+1;
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end

%% Parallel loop
parfor k=1:length(position)
    
    getFRETDataIX(position{k},bgdir,rawdir,datadir);

end
disp('done!');

%% Bleaching correction
% calculate average bleaching from all movies in the same 
% datadir='K:\data\150617_NG\rawdata\RhoA2G_20s_2\data_150618';
close;
position={}
k=0;
for row=2:3
    for col=2
        for site=1:3
            k=k+1;
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end


bleach_raw_all=[];
for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
%     bleach_raw(bleach_raw>3)=NaN;
%     bleach_raw(bleach_raw<1)=NaN;
    bleach_raw_all(:,k)=bleach_raw/nanmedian(bleach_raw);
    %bleach_raw_all=bleach_raw/median(bleach_raw);
end
bleach_mean=nanmedian(bleach_raw_all,2);
plot(1:length(bleach_mean),bleach_raw_all);hold on
plot(1:length(bleach_mean),bleach_mean);

timepts=1:length(bleach_mean);
[fitvals]=polyfit(timepts',bleach_mean,1);
corr=polyval(fitvals,timepts);

plot(timepts,polyval(fitvals,timepts),'LineWidth',3);axis([0 length(bleach_mean) 0 1.5]);

%% correct bleaching
slope=fitvals(1);
slope=0;

close;
position={}
k=0;
for row=2:3
    for col=2
        for site=1:3
            k=k+1;
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end



k=1:length(position)
for k=1:length(position) 
    
%     load([datadir,filesep,position{k},'_Bleach_raw.mat']);
%     bleach_mean=bleach_raw;
%     timepts=1:length(bleach_mean);
%     [fitvals]=polyfit(timepts,bleach_mean,1);
%     corr=polyval(fitvals,timepts);  
%     slope=fitvals(1);
%         
    correctBleachingLinear(position{k},slope,datadir);

end
disp('done!');
    

