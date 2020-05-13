%% ParallelFRET.m
% This script sets up parallel processing to obtain FRET data of time-lapse
% series by calling getFRETData in a parfor loop.
% Arnold, 23 May 2015

%% Parameter setup
clear 
clc 
close all
%warning off
root = 'K:\data\150424_IX\RhoA2G-HT73-30s';
rawdir=[root,filesep,'rawdata'];
bgdir='K:\data\150424_IX\background-bin1';
datadir=[root,filesep,'data_150527'];
mkdir(datadir);

rows=3;
cols=2:6;   %[3 7]
sites=1:2; %1:9

% %%
% shot={};
% i=0;
% for rows=4
%     for cols=2:3
%         for sites=1:2
%             i=i+1
%             shot{i}=[num2str(rows),' ',num2str(cols),' ',num2str(sites)];
%             shot=shot';
%         end
%     end 
% end
%%
manualwells = [
    1 1 1;
    
    ];

manualcontrol=0;
numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
shots=numrows*numcols*numsites;


if manualcontrol==1
    shots=size(manualwells,1);
end
time1=tic;
parfor shot=1:shots
    if manualcontrol==1
        row=manualwells(shot,1);
        col=manualwells(shot,2);
        site=manualwells(shot,3);
    else
        siteidx=mod(shot,numsites);
        if siteidx==0
            siteidx=numsites;
        end
        site=sites(siteidx);
        colidx=mod(ceil(shot/numsites),numcols);
        if colidx==0
            colidx=numcols;
        end
        col=cols(colidx);
        rowidx=ceil(shot/(numcols*numsites));
        row=rows(rowidx);
    end
    %fprintf([num2str(row),'_',num2str(col),'_',num2str(site),'\n']);
    
    %%% Format Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FormatFiles(row,col,site);
    
    %%% Timelapse %%%%%%%%%%%%%%%%%%%%%%%%%
    getFRETDataIX(row,col,site);
    
    %%% Immunostain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Immunostain_1_CalcBleedthroughRate(row,col,site);
    %Immunostain(row,col,site);
    %Immunostain_2_AddToTimelapse(row,col,site);
end

toc(time1)