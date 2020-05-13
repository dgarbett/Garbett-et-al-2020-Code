% combine_pre_post_150506.m
% combine files of images acquired pre and post drug addition, give them
% one consecutive index and move them to a new folder. 
clear; clc;
predir='K:\data\150701_IX\pre3';
postdir='K:\data\150701_IX\post3';
targetdir='K:\data\150701_IX\CR_thrombin';
if ~exist(targetdir)
    mkdir(targetdir);
end
preframes=10;
postframes=60;

%% move predrug files
for row=2:3
    for col=8:9
        for site=4
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            for frame=1:preframes
                sourcefile_w1=[predir,filesep,shot,'_CFP_',num2str(frame),'.tif'];
                sourcefile_w2=[predir,filesep,shot,'_FRET_',num2str(frame),'.tif'];
                targetfile_w1=[targetdir,filesep,shot,'_CFP_',num2str(frame),'.tif'];
                targetfile_w2=[targetdir,filesep,shot,'_FRET_',num2str(frame),'.tif'];
                movefile(sourcefile_w1,targetfile_w1,'f');
                movefile(sourcefile_w2,targetfile_w2,'f');
                disp([shot,'_',num2str(frame)]);
            end
        end
    end
end

%% move postdrug files
for row=2:3
    for col=8:9
        for site=4
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            for frame=1:postframes
                fileind=frame+preframes;
                sourcefile_w1=[postdir,filesep,shot,'_CFP_',num2str(frame),'.tif'];
                sourcefile_w2=[postdir,filesep,shot,'_FRET_',num2str(frame),'.tif'];
                targetfile_w1=[targetdir,filesep,shot,'_CFP_',num2str(fileind),'.tif'];
                targetfile_w2=[targetdir,filesep,shot,'_FRET_',num2str(fileind),'.tif'];
                movefile(sourcefile_w1,targetfile_w1,'f');
                movefile(sourcefile_w2,targetfile_w2,'f');
               disp([shot,'_',num2str(frame)]);
            end
        end
    end
end