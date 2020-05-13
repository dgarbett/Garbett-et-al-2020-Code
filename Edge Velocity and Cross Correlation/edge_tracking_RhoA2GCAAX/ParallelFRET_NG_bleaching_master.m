%% ParallelFRET_NG.m
% This script sets up parallel processing to obtain FRET data of time-lapse
% series by calling getFRETData in a parfor loop.
% For all movies in the same folder, a 
% Arnold, 9 June 2015

%% Parameter setup
clear;clc;close all; 
root = '/Users/abisaria/Desktop/20151113_(HeLaTiamFRBRap)';
rawdir=[root,filesep,'Rapamycin_B8_test_1'];
bgdir=[root,filesep,'background'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[rawdir,filesep,'data'];
if ~exist(datadir)
    mkdir(datadir);
end

jitter=0; %jitter correction, 0 or 1.
startpos=0;
endpos=5;

j=0;
for i=startpos:endpos
    j=j+1;
    position{j}=['Pos',num2str(i)];
end

%position={'Pos1' 'Pos3' 'Pos4' 'Pos5'}

%% loop
for k=1:length(position)
    
    %getFRETDataNG(position{k},bgdir,rawdir,datadir);
    interactive_edge_tracking(position{k},bgdir,rawdir,datadir,jitter);
end
disp('done!');

%% Bleaching correction
% calculate average bleaching from all movies in the same 
% datadir='K:\data\150617_NG\rawdata\RhoA2G_20s_2\data_150618';
close;
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
%[fitvals]=polyfit(timepts',bleach_mean,1);
fitpara=fit(timepts',bleach_raw_all(:,1),'exp2');
%plot(f,timepts',bleach_raw_all(:,3));axis([0 length(bleach_mean) 0.8 1.2]);

% [fitvals]=polyfit(timepts',bleach_raw_all(:,3),1);
% corr=polyval(fitvals,timepts);

%plot(timepts,polyval(fitvals,timepts),'LineWidth',3);axis([0 length(bleach_mean) 0.8 1.2]);



%% correct bleaching
%slope=fitvals(1);
%parfor k=1:length(position)
for k=1:length(position) 
    
%     load([datadir,filesep,position{k},'_Bleach_raw.mat']);
%     bleach_mean=bleach_raw;
%     timepts=1:length(bleach_mean);
%     [fitvals]=polyfit(timepts,bleach_mean,1);
%     corr=polyval(fitvals,timepts);  
%     slope=fitvals(1);
%         
    correctBleachingExp2(position{k},fitpara,datadir);

end
disp('done!');
    

