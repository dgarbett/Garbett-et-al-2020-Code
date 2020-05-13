%% Test protrusion/FRET measurements
clear; close all; clc;

%% Parameters for cell edge parametrization
nFretWindows=360;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=5;                                      % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge
edgeDepthDist=15;                                   % Number of pixels deep for the windows for computing FRET values.
%dataFolder=[root filesep 'Data'];
startFrame=1;
endFrame=60;
binning=1;

%colormap('parula');

%% Load previously determined sequence of masks
%load('D:\Stanford\Matlab\XCorr\From Arnold For cross-correlation between FRET and edge dynamics (150109)\FRET_analysis\sequence\masks and coors.mat');
load('K:\data\150423_NG\data_150527\RhoA2G_1_RatioData.mat');

% load('K:\data\141121_NG\data_150224\CC-Cdc42_15s_5_pos_0_RatioData.mat');
% perform tracking of mask centroid
if min(cellfun(@(x) size(x,1),cellCoors))>0
    traj=ultTrackAnnSearch(cellCoors,'pairrule','fwdbckmtch','maxdisp',100,'verbose',false);
else
    traj={};
end
fprintf('%i trajectories.\n',length(traj));


%% Display mask and number objects
imagesc(maskFinal{1});hold on;
for n=1:length(traj);
    text(double(traj{n}(1,1)),double(traj{n}(1,2)),num2str(traj{n}(1,4)));
end
 
%%
thisTraj=traj{17};

frameCount=0;
 for imnum=startFrame:endFrame
            frameCount=frameCount+1;
            
            % Make a one cell mask
            objects=regionprops(maskFinal{imnum},'PixelIdxList','PixelList','Centroid','BoundingBox');
            cellCent=round(thisTraj(find(thisTraj(:,end)==imnum,1),1:2));
            cellNum=find(round(arrayfun(@(x) x.Centroid(1),objects))==cellCent(1) & round(arrayfun(@(x) x.Centroid(2),objects))==cellCent(2));
            thisMask=false(size(maskFinal{imnum}));
            thisMask(objects(cellNum).PixelIdxList)=true;        
            thisMask=false(size(maskFinal{imnum})); % sets a frame of the size of the mask to zero
            thisMask(objects(cellNum).PixelIdxList)=true; % sets values of object 2 to 1
                       
            % Parametrize cell edge and compute protrusion values
            if frameCount==1
                [edgeCoors{frameCount} edgeCoorsSmoothed{frameCount}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning));
            else
                [edgeCoors{frameCount},edgeCoorsSmoothed{frameCount}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning),edgeCoors{frameCount-1});
                windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                protvals(:,frameCount-1)=vect(computeProtrusionValues(edgeCoors{frameCount-1},edgeCoorsSmoothed{frameCount-1},edgeCoors{frameCount}));
            end
            windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
            
            % Compute edge regions and local fret ratio values
            labelMask=getWindowLabelMap(thisMask,windowCoors{frameCount},edgeDepthDist/binning);
            for k=1:size(windowCoors{frameCount},1)
                %fretvals(k,frameCount)=sum(imFRET(labelMask==k))/sum(imCFP(labelMask==k));
                fretvals(k,frameCount)=mean(imRatio{frameCount}(labelMask==k));
            end
            disp(num2str(imnum));
        end
        %% Visualize edge coordinates
%         close all;
%         plot(edgeCoors{1}(:,1),edgeCoors{1}(:,2),'-r');hold on
%         plot(edgeCoorsSmoothed{1}(:,1),edgeCoorsSmoothed{1}(:,2));
%         %%
        protvalsWindow=zeros(size(fretvals)-[0 1]);
        for k=1:edgeOversamplingParam
            protvalsWindow=protvalsWindow+protvals(k:edgeOversamplingParam:end,:);
        end
        % Normalize protrusion data by dividing by the edgeOversampling factor and dividing by the time between images.
        %protvalsWindow=(protvalsWindow/edgeOversamplingParam)./repmat(diff(timeData(startFrame:endFrame))',[nFretWindows 1]);

protRange=[round(prctile(protvalsWindow(:),2),1),round(prctile(protvalsWindow(:),98),1)];        
subplot(1,2,1);imagesc(protvalsWindow,protRange);
subplot(1,2,2);imagesc(fretvals);


disp('done!');