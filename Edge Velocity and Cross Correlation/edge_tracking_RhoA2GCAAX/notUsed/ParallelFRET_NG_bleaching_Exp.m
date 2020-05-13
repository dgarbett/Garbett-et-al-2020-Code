%% ParallelFRET_NG.m
% This script sets up parallel processing to obtain FRET data of time-lapse
% series by calling getFRETData in a parfor loop.
% For all movies in the same folder, a 
% Arnold, 9 June 2015

%% Parameter setup
clear;clc;close all; 
root = 'K:\data\150620_NG\CDH5-Ruby_CC-Rac1_1min_4';
rawdir=[root];
bgdir=['K:\data\150620_NG',filesep,'background'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[rawdir,filesep,'data_150622'];
if ~exist(datadir)
    mkdir(datadir);
end

startpos=0;
endpos=5;

j=0;
for i=startpos:endpos
    j=j+1;
    position{j}=['Pos',num2str(i)];
end

%% Parallel loop
parfor k=1:length(position)
    
    getFRETDataNG_RFP(position{k},bgdir,rawdir,datadir);

end
disp('done!');

%% Bleaching correction
% calculate average bleaching from all movies in the same 
% datadir='K:\data\150617_NG\rawdata\RhoA2G_20s_2\data_150618';
close;
for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
%     bleach_raw(bleach_raw>3)=NaN;
%     bleach_raw(bleach_raw<1)=NaN;
    bleach_raw_all(:,k)=bleach_raw;
end
bleach_mean=nanmedian(bleach_raw_all,2);
plot(1:length(bleach_mean),bleach_raw_all);hold on
plot(1:length(bleach_mean),bleach_mean);

timepts=1:length(bleach_mean);

% expfit=fit(timepts',bleach_mean,'exp1');
% plot(expfit,timepts,bleach_mean);

[fitvals]=polyfit(timepts',bleach_mean,1);
corr=polyval(fitvals,timepts);

plot(timepts,polyval(fitvals,timepts),'LineWidth',3);axis([0 length(bleach_mean) 1 3]);

%% correct bleaching
slope=fitvals(1);
%slope=0;
%parfor k=1:length(position)
for k=1:length(position) 
    
%     load([datadir,filesep,position{k},'_Bleach_raw.mat']);
%     bleach_mean=bleach_raw;
%     timepts=1:length(bleach_mean);
%     [fitvals]=polyfit(timepts,bleach_mean,1);
%     corr=polyval(fitvals,timepts);  
%     slope=fitvals(1);
%         
    correctBleaching(position{k},slope,datadir);

end
disp('done!');
    

