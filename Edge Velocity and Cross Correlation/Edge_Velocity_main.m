%Image processing for/ratio imaging
%
% Input
%   - Movie sequence with one image per timepoint, acquired using
%   beamsplitter
%   - A background image acquired the same way, ideally an average of
%   multiple background images.
%   - Images of the alignment grid to calculate the offset parameters
%
% Output
%   - movie file with ratio images
%   - cell array with ratio images stored
%
% Background subtraction (subBG.m)
%   - After generating a mask containing FRET-probe expressing cells, the
%   background level is determined from background regions and normalized
%   by using the background image to account for uneven illumination
%   (flatfield correction)
%
% Arnold Hayer, 150115
%
% Adapted functions for use of the Hamamatsu sCMOS camera.


%% Working folder
clear
clc
close all
warning off

%what direction do file separations use
filesep = '/';

%Files
root = '/Volumes/AB_Data_2/DG_analysis/Edge Velocity/Stuff for Edge Velocity';
rawdir=[root,filesep];
datadir=[root,filesep,'data'];

prefix = 'Capture'; %prefix for loading files based on condition or prefix from 3i
%%RENAME FOR NEW ANALYSIS <-DAMIEN
output = 'tracking_20181002_30windows';%  
moviebin = 2; %binning (1 or 2)
umPerPx = .325;
visualOutput = 0; %output some intermediate movies?
RGBwrite = 1; %needed for tracking code, keep at 1

runTracking = 0; %run cell tracking
RunGradientAnalysis=0;
RunFullGradientAnalysis = 0;
RunEdgeVelocity = 1;

%Edge Velocity Analysis
% diff number of windows <-DAMIEN
nFretWindows=80;                                 % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=10; %window size for edge

%channel names (can be C1/or actual name - note C0 is the membrane for
%averaging, al are used for segmentation in this case).


channelNames = {'FT' 'PLS3'};
channelCombs = [2 1];
 
comboNames = {'Norm PLS3'};
plotIndexes=[1 2];



timeStep = 2; %minutes between time steps;
window = 3; %specify pixel window for how much of a ring you want if at all.
%if you want to be the same as edge velocity analysis make same
%as edgeDepthDist
numhist = 20;%how many bins to calculate gradient across
allHist = numhist*3;

allNames = [channelNames comboNames];
segmentationCh = [1];
if length(comboNames) > size(channelCombs,1)
    error('Need extra channel combo');
elseif size(channelCombs,1) > length(comboNames)
    for i = length(comboNames)+1: size(channelCombs,1)
        comboNames{i} = 'N/A';
    end
end


RGBdir = ([root,filesep,output,filesep,'RGB',filesep]);
if ~exist(RGBdir,'dir')
    mkdir(RGBdir);
end
if ~exist(datadir,'dir')
    mkdir(datadir);
end
cd(rawdir)
outputdir = ([root,filesep,output,filesep]);
if ~exist(outputdir,'dir')
    mkdir(outputdir);
end



%%% set up other variables, based on 40x.  change accordingly
    ring = ceil(window/umPerPx);
    nucr = ceil(10/umPerPx);
    debrisarea = ceil(400/(umPerPx^2));
    boulderarea = ceil(2100/(umPerPx^2));
    memdebrisarea=debrisarea;
    coff = ceil(50/umPerPx);


for c = 1:length(channelNames)
    bg_raw(:,:,c) = double(imread(['background',filesep,'AVG_bgC',num2str(c-1),'.tif']));
end
%%
for ii = 1:length(channelNames)
    files = strcat('*',prefix, '*C', num2str(ii-1),'*.tif');
    IMChStack{ii} = dir(files);
    
end
numFiles = numel(IMChStack{ii});


%%
for I = [1];
    allMem = [];
    for ii = 1:length(channelNames)
        
        IMChStackI = IMChStack{1,ii};
        IMChName{ii} = IMChStackI(I).name;
        IMInfo{ii} = imfinfo(IMChStackI(I).name);
        
        filename = regexprep(IMChName{ii},['_C', num2str(ii-1),'.tif'],'');
        thisStack = IMChStack{1,segmentationCh(1)};
        thisInfo = IMInfo{ii};
        chMem = [];
        for f = 1:numel(thisInfo)
            im_mem = double(imread(IMChName{ii}, f, 'Info', thisInfo));
            chMask = segment_logMultiThresh(im_mem,debrisarea);
            chMask =  imdilate(chMask,strel('disk',10));
            im_mem(chMask ==1) = nan;
            chMem(:,:,f) = im_mem;
        end
        meanCh = nanmean(chMem,3);
        PixMinMax=double([round(prctile(meanCh(:),1)) round(prctile(meanCh(:),99))]);
        IntStep=((PixMinMax(2)-PixMinMax(1))/100);  
        [f,xi]=ksdensity(meanCh(:),PixMinMax(1):IntStep:PixMinMax(2));
        [M, mi] = max(f);
        bg_value(ii) = xi(mi);  
    end
    

    startFrame = 1;%assume start frame is 1 unless other specified
    endFrame = numel(IMInfo{1,ii}); % image until end of stack
    
    frames = startFrame:4:endFrame;
    badframes = NaN(endFrame,1);
    
    %tracking parameters
    jitters = zeros(endFrame,2);
    jitterwindow=3*nucr;
    blobthreshold = -0.02;  %default 10xbin1=-0.02 20xbin2=-0.03;
    memblocksize = 500;
    maxmemnum = 100;
    maxjumpmem = nucr;
    %1:x_loc 2:y_loc 3.size? 4.? 5? 6? 7? add in nuc x,y and other things
    
    
    
    if runTracking == 1
        %set up trace data files
        memtracedata = NaN(maxmemnum,endFrame,5);
        memtracking = NaN(maxmemnum,5);
        
        %set up initial tracking parameters form first 3 images
        for f = 1:3
            thisStack = IMChStack{1,segmentationCh(1)};
            thisInfo = IMInfo{segmentationCh(1)};
            im_mem = double(imread(IMChName{segmentationCh(1)}, f, 'Info', thisInfo));
            
            if length(segmentationCh) > 1
                for ii = 1:length(segmentationCh)-1
                    thisInfo = IMInfo{segmentationCh(ii+1)};
                    nextCh = double(imread(IMChName{segmentationCh(ii+1)}, f, 'Info', thisInfo));
                    im_mem = im_mem+nextCh;
                    mem_temp(:,:,f) = im_mem;
                end
            else
                mem_temp(:,:,f) = double(imread(IMChName{segmentationCh}, f, 'Info', thisInfo));
            end
        end
        
        height = thisInfo(1).Height;
        width = thisInfo(1).Width;
        %tracking initialization
        [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width] = timelapsesetup_membrane(mem_temp,nucr,debrisarea,badframes,RGBwrite);
        
        
        folder = [rawdir,filesep];
        
        %go through all frames
        for i = length(frames)%firstgoodindex:endFrame
            f = frames(i);
            
            %%PROBABLY CHANGE THIS
            %reads in each tiff stack
            
            for ii = 1:length(channelNames)
                thisInfo = IMInfo{ii};
                allChRaw(:,:,ii) = double(imread(IMChName{ii}, f, 'Info', thisInfo));
            end
            %allChRaw =
            rawMem = [];
            rawMem = allChRaw(:,:,segmentationCh(1));
            %tracking based on 2 channels
            %for ii = 1:length(segmentationCh)-1
               % rawMem = rawMem + allChRaw(:,:,segmentationCh(ii+1));
           % end
            Ch0_raw = allChRaw(:,:,1);
            Ch1_raw = allChRaw(:,:,2);
            bgmask=getBGMask(Ch0_raw+ Ch1_raw);
            
            
            imCh0bg=subBG(Ch0_raw,bgmask,bg_raw(:,:,1));
            imCh1bg=subBG(Ch1_raw,bgmask,bg_raw(:,:,2));
           
           
            sumImage=Ch0_raw+Ch1_raw;
            mem_mask=getCellMask_forDG(sumImage,4000);
            %em_mask = 
            [mem_label,nummem] =  bwlabel(mem_mask); %number cell objects
%             if sum(mem_label(:)) == 0 %if nothing was tracked, note it as a bad frame, and write empty mask
%                 badframes(f) = 1;
%                 if f == 0
%                 imwrite(uint16(mem_label),[RGBdir,filename,'_memedge.tif']);
%                 else                  
%                 imwrite(uint16(mem_label),[RGBdir,filename,'_memedge.tif'],'WriteMode','append');
%                 end
%                 continue
%             end
            
            %get info about centroid/cell intensity
            mem_info = struct2cell(regionprops(mem_mask,rawMem,'Area','Centroid','MeanIntensity')');
            mem_area = squeeze(cell2mat(mem_info(1,1,:)));
            mem_center = squeeze(cell2mat(mem_info(2,1,:)))';
            mem_density = squeeze(cell2mat(mem_info(3,1,:)));
            mem_mass = mem_density.*mem_area;
            
            %put information into memcurdata
            if numel(mem_area) == 0
            elseif numel(mem_area)== 1
                memcurdata = [mem_center(1),mem_center(2),mem_area,mem_mass'];
            else
                memcurdata = [mem_center(:,1),mem_center(:,2),mem_area,mem_mass];
            end
            
            
            if i>1
                lastgoodframe = find(badframes==0,1,'last'); %infd the last "good frame"
                H2BvsNLS = 1; %1:H2B 2:NLS ... not important for us
                
                %track based on last good frame
                [memtracedata,memcurdata,memtracking,mem_label] = adaptivetrack_membrane_mass(lastgoodframe,f,memtracedata,memcurdata,memtracking,rawMem,mem_label,nucr,jitters(f,:),maxjumpmem,H2BvsNLS);
                badframes(f)=0;
            end
            
            %write the masks
            extractmemmask = bwmorph(mem_label,'remove');
            %mem_label_ring = mem_label - imerode(mem_label,strel('disk',1,0));
            %write output files after bg subtraction
            
            if f == 1
                imwrite(uint16(mem_label),[RGBdir,filename,'_memedge.tif']);
            else
                imwrite(uint16(mem_label),[RGBdir,filename,'_memedge.tif'],'WriteMode','append');
            end
            
            
            %cell and membrane tracking
            memid = find(~isnan(memcurdata(:,1)));
            numtrackmem = numel(memid);
            
            %get centroid and area information for membranes
            memcurdata = memcurdata(memid,:);
            mem_center = memcurdata(:,[1 2]);
            mem_area = memcurdata(:,3);
            mem_mass = memcurdata(:,4);
            mem_info = regionprops(mem_label,'PixelIdxList');
            %mem_raw1 = regionprops(mem_label,imMem_raw,'PixelValues');
            mem_raw2 = regionprops(mem_label,rawMem,'PixelValues');
            mem_nuc_num =[];
            
            memIDvec = {mem_info.PixelIdxList};
            memidvec = ~cellfun(@isempty,memIDvec );
            
            mem_info = mem_info(memidvec);
            %mem_raw1 = mem_raw1(memidvec);
            mem_raw2 = mem_raw2(memidvec);
            %need to make sure only iterate through mem_info if it isn't
            %empty
            if size(mem_center,1) == 1 && ~isempty(mem_raw2)
                mem_nuc_num(1,1) = mode(mem_raw2.PixelValues(mem_raw2.PixelValues>0));
            elseif ~isempty(mem_raw2)
                for m = 1:size(mem_center,1)
                    mem_nuc_num(m,1) = mode(mem_raw2(m).PixelValues(mem_raw2(m).PixelValues>0));
                end
            end;
            
            %add back to memtrace data
            memtracedata(memid,f,:) = [mem_center,mem_area,mem_mass,mem_nuc_num];
            if maxmemnum - max(memid) < memblocksize
                memtempdata = NaN(memblocksize,endFrame,5);
                memtemptrack = NaN(memblocksize,5);
                memtracedata = [memtracedata;memtempdata];
                memtracking = [memtracking;memtemptrack];
                maxmemnum = maxmemnum+memblocksize;
            end
        end;
        
        %         if sum(badframes)+2 >= size(memtracedata,2)
        %             save([outputdir,'memtracedata_',filename,'.mat'],'memtracedata');
        %             continue;
        %         end;
        
        %link cell tracks
        [memtracedata,~,~] = postprocessing_nolinking(memtracedata,memid,jitters,badframes,memtracking,maxmemnum,nucr);
        meminfo = imfinfo([RGBdir,filename,'_memedge.tif']);
        
        %%remove all NaN rows
        for ii = 1:size(memtracedata,2)
            
            memEdgeImage=imread([RGBdir,filename,'_memedge.tif'],(ii/4)+1, 'Info', meminfo);
            cellsinFrame = unique(memEdgeImage(memEdgeImage>0));
            if isempty(cellsinFrame)
                memtracedata(:,ii,:) = NaN;%*ones(size(memtracedata,1),1,size(memtracedata,3));
            else
                cellsVerified = ismember((1:size(memtracedata,1)),cellsinFrame);
                memtracedata(~cellsVerified,ii,:) = NaN;
            end
            
        end
        
        save([outputdir,'memtracedata_',filename,'.mat'],'memtracedata');
        
    end %runTracking ==1;
    
    
    
    % start ratio analysis
    close all;
    
    if runTracking == 0 %if you don't rerun tracking load the tracedata file
        load([outputdir,'memtracedata_',filename,'.mat']);
    end
    
    %run through tracked cells
    for m = 1:size(memtracedata,1)
        saveFile = 1; %switch in case a cell sucks and we don't want to output the data
        
        clear edgeCoors protvals windowCoors shape edgevals;
        edgevals =[];
        plotTitle = 'Normalized Intensity'; %what to label plot y-axis
        
        
        meminfo = imfinfo([RGBdir,filename,'_memedge.tif']);% load tracking tiff
        cellData=[];
        
        
        cellData = squeeze(memtracedata(m,:,:));
        disp(['cell = ',num2str(m)]);
        
        %see that the cell was tracked
        longestrun = regionprops(~isnan(cellData(:,1)),'Area','PixelIdxList');
        if isempty(longestrun)
            continue;
        end
        
        %find the longest tracked run for the biggest object found for a cell (in case
        %something small is left over after erosion etc)
        [sortedAreas, sortIndexes] = sort([longestrun.Area], 'Descend');
        idx_first = longestrun(sortIndexes(1)).PixelIdxList(1);
        idx_last = longestrun(sortIndexes(1)).PixelIdxList(end);
        
        
        %make sure nothign is funky with tracking and that it was tracked
        %for at least 5 frames
        if idx_last > numel(meminfo)
            idx_last = numel(meminfo);
        end
        if idx_last-idx_first < 5
            continue;
        end
        
        frames_tracked = idx_first:idx_last;
        %initalize matrices;
        edgeKymograph = [];
        distKymograph=[];
        pxIntRatio=[];
        
        %read in the 1st tiff
        memEdgeImage=imread([RGBdir,filename,'_memedge.tif'],idx_first, 'Info', meminfo);
        memEdgeImageCell=memEdgeImage==m ;

            %calculate the crop area
            %change for subamplings <-DAMIEN
            subsample = 1:2:length(frames_tracked)-1;
            for frames=1:length(subsample) %iterate through tracked frames;
                frame = subsample(frames);
                %read in membrane edge for image;
                memEdgeImage=imread([RGBdir,filename,'_memedge.tif'],frame, 'Info', meminfo);
                memEdgeImageCell=memEdgeImage==m ;
                memEdgeImageCell = imerode(memEdgeImageCell,strel('disk',5));
                memEdgeFill =memEdgeImage==m ;
                
                %make a ring of size "ring" in this case 3 um
                cellInside = double(imerode(memEdgeImageCell,strel('disk',ring)));
                memEdgeImageCell = logical(memEdgeImageCell.* ~cellInside);
                memEdgeImageCell = bwareaopen(memEdgeImageCell, 20);
                %memEdgeImageCell = imopen(memEdgeImageCell, strel('disk',1));
                memEdge=bwperim(memEdgeImageCell);
                
                
                for ii = 1:length(channelNames)
                    thisInfo = IMInfo{ii};
                    chRaw = double(imread(IMChName{ii}, frame, 'Info', thisInfo));
                    %ChBG = nanmean(vect(double(immultiply(chRaw,~memEdgeFill))));
                    %ChBG = nanmin(chRaw);
                    chRaw = chRaw - bg_raw(ii);
                    %rangeVals = range(chRaw(:));
                    %chRaw = rescale(chRaw,1,rangeVals);
                    %chRaw = imtophat(chRaw,strel('disk',400));
                    chRaw(memEdgeFill ~= 1)  = nan;    
                    chRaw=ndnanfilter(chRaw,fspecial('disk',2),'replicate'); %changed disk to 2 to blur less
                   % medNormChannel = (chRaw./nanmedian(chRaw(:)));
%                     if nanmin(medNormChannel(:)) <= 0
%                         medNormChannel = medNormChannel - (nanmin(medNormChannel(:))*1.00001);
%                     end
                     allChNorm(:,:,ii) =  chRaw;
                    
                end
                
                for ii = 1:size(channelCombs,1)
                    jj = ii+length(channelNames);
                    ind1 = channelCombs(ii,1);
                    ind2 = channelCombs(ii,2);
                    ch1 = ndnanfilter(allChNorm(:,:,ind1),fspecial('average',3),'replicate');
                    ch2 = ndnanfilter(allChNorm(:,:,ind2),fspecial('average',3),'replicate');
                    
                    allChNorm(:,:,jj) = ch1./ch2;
                 
                end

                if sum(memEdge(:)) > 0
                    %start analzying data
                    stats = regionprops(memEdgeFill, 'Centroid');
                    memOutline = regionprops(memEdge,'PixelList');
                    sizeVals = regionprops(memEdgeFill,'Eccentricity','Area','PixelList');
                    
                    PxLocs = regionprops(memEdgeImageCell,'PixelList','PixelIdxList');
                    
                    allInts = [];
                    for kk = 1:size(allChNorm,3)
                        chVals = regionprops(memEdgeImageCell,allChNorm(:,:,kk),'PixelValues');
                        pxVals = chVals.PixelValues;
                        allInts = [allInts pxVals];
                    end
                    
                    centroid = stats.Centroid;
                    RingPx = PxLocs.PixelList;
                    AllPx = sizeVals.PixelList;
                    shape(1:2,frames) = [sizeVals.Eccentricity sizeVals.Area];
                    
                    
                    if RunEdgeVelocity
                        %%Start Edge Velocity Analysis
                        cellMaskVals = regionprops(memEdgeFill,'PixelIdxList');
                        
                        if numel(cellMaskVals) > 1
                            maskEl = 0;
                            lengthEl = 0;
                            for ii = 1:numel(cellMaskVals)
                                if length(cellMaskVals(1).PixelIdxList) > lengthEl;
                                    maskEl = ii;
                                    lengthEl = length(cellMaskVals(1).PixelIdxList);
                                end
                            end
                            
                            cellMaskVals = cellMaskVals(maskEl);
                        end
                        cellCent=round(centroid);
                        thisMask=false(size(memEdgeFill));
                        thisMask(cellMaskVals.PixelIdxList)=true;
                        
                        % Parametrize cell edge and compute protrusion values
                        if frames==1
                            [edgeCoors{frames} edgeCoorsSmoothed{frames}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/moviebin));
                            protvals=[];
                        else
                            [edgeCoors{frames} edgeCoorsSmoothed{frames}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/moviebin),edgeCoors{frames-1});
                            %windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                            protvals(:,frames-1)=vect(computeProtrusionValues(edgeCoors{frames-1},edgeCoorsSmoothed{frames-1},edgeCoors{frames}));
                        end
                        windowCoors{frames}=edgeCoors{frames}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                        
                        % Compute edge regions and local fret ratio values
                        labelMask{frames}=getWindowLabelMap(thisMask,windowCoors{frames},ring);
                        
                        for k=1:size(windowCoors{frames},1)
                            %fretvals(k,frames)=sum(imFRET(labelMask==k))/sum(imCFP(labelMask==k));
                            for jj = 1:size(allChNorm,3)
                                thisIM = allChNorm(:,:,jj);
                                edgevals(k,frames,jj) =nanmean(thisIM(labelMask{frames}==k));
                            end
                        end %edgeVals for each channel
                    end %if run edge velocity
                else
                    saveFile = 0;
                    break;
                end
            end % through tracked frames
            
            % smooth prot values
            if RunEdgeVelocity
                if length(protvals) > 0
                    protvalsWindow=zeros([size(edgevals,1) size(edgevals,2)]-[0 1]);
                    for k=1:edgeOversamplingParam
                        protvalsWindow=protvalsWindow+protvals(k:edgeOversamplingParam:end,:);
                    end
                end
                
                protvalsWindowF=ndnanfilter(protvalsWindow,fspecial('disk',2),'replicate');
                
                edgeSensorVals = NaN(size(protvalsWindowF,1),size(protvalsWindowF,2)+1,size(allChNorm,3));
                for jj = 1:size(allChNorm,3)
                    temp = squeeze(edgevals(:,:,jj));
                    infvals = temp ==Inf;
                    negvals = temp <= 0;
                    
                    temp(infvals==1) = NaN;
                    temp(negvals == 1) = NaN;
                    edgeSensorVals(:,:,jj) = ndnanfilter(temp,fspecial('disk',2),'replicate');
                end
            end
            if saveFile
                if RunEdgeVelocity && RunGradientAnalysis
                    save([output,'/',filename,'_cell_analysis_v2',num2str(m)],'MeanHistTime','MedianHistTime','cellData','v_all','a_all' ,'protvals','protvalsWindow','protvalsWindowF','edgeSensorVals','MSD','shape','edgevals','labelMask');
                elseif RunGradientAnalysis
                    save([output,'/',filename,'_cell_analysis_v2',num2str(m)],'MeanHistTime','MedianHistTime','cellData','v_all','a_all' ,'MSD','shape');
                elseif RunEdgeVelocity
                    save([output,'/',filename,'_cell_analysis_v2',num2str(m),'.mat'],'cellData','protvals','protvalsWindow','protvalsWindowF','edgeSensorVals','shape','labelMask');
                else
                    save([output,'/',filename,'_cell_analysis_v2',num2str(m)],'cellData','v_all','a_all','MSD','shape');
                end
            end
            close all;
   
        
        disp(['End of cell ',num2str(m)]);
    end%through tracked membranes
end  %through imagges

%%


%figure out where the cell is moving?
timeStep = 1;
longestrun = regionprops(~isnan(cellData(:,1)),'Area','PixelIdxList');

%find the longest tracked run for the biggest object found for a cell (in case
%something small is left over after erosion etc)
[sortedAreas, sortIndexes] = sort([longestrun.Area], 'Descend');
idx_first = longestrun(sortIndexes(1)).PixelIdxList(1);
idx_last = longestrun(sortIndexes(1)).PixelIdxList(end);
idx_first = 1;

a_all = NaN(1,idx_last-idx_first);
frames_tracked = idx_first:idx_last;

%calculate a "fitted" x,y position basde on fitted trajectory
dAll = [];
%calculate the speed and trajectory
for frames=1:length(frames_tracked)-2
    frame = frames_tracked(frames);
    
    %speed is valcualted based on actual speed
    tm=[cellData(frame,1) cellData(frame,2)];
    tm_step=[cellData(frame+2,1) cellData(frame+2,2)];
    d=tm_step-tm; % Vectors connecting two successive points.
    [atot, vtot]=cart2pol(d(:,1),d(:,2));
    dAll(:,frames) = d; 
    v_all(:,frames)=vtot;
    
    %trajectory is based on fitted tracking
    
    [atot, vtot]=cart2pol(d(:,1),d(:,2));
    a_all(:,frames)=rad2deg(atot);
end %* end angle/velocity calculation;


timestart = 1;
timeend = size(MeanHistTime,2);

xlocs2 = ([1:10:size(MeanHistTime,2)]-1); 
 xlabels2 = arrayfun(@(i) num2str(ceil(i*timeStep)),xlocs2,'UniformOutput',false);

a =  makeColorMap([25 025 112]/255,[1 1 1 ], [139 0 0 ]/255,600);
%%
figure;


ax2 =  subplot(1,3,1);
imagesc(imgaussfilt(protvalsWindowF*umPerPx/timeStep,1));

set(gca, 'XTick',xlocs2, 'Xticklabel',xlabels2);
title('Protrusion/Retraction  (um/hr)');
ylabel('Position');
colormap(ax2,a)
maxvel = prctile(abs(vect(protvalsWindowF*umPerPx/timeStep)),92); 
caxis(ax2,[-maxvel  maxvel ]);
colorbar('EastOutside');
xlabel('Time (min)');
xlim([timestart timeend-1]);
set(gca, 'XTick',xlocs2, 'Xticklabel',xlabels2);
 
 
 
ax3=subplot(1,3,2);
cact = edgeSensorVals(:,:,1);
imagesc(cact)%./edgeSensorVals(:,:,4));
colormap(ax3,parula);
climits = [prctile(vect(cact),5) prctile(vect(cact),98)];
%climits=[.5 5];
caxis(ax3,climits);
title('F-tractin');
 set(gca, 'XTick',xlocs2, 'Xticklabel',xlabels2);
 set(gca, 'YTick',[], 'Yticklabel',[]);
colorbar('EastOutside');
set(findall(gcf,'-property','FontSize'),'FontSize',12);
xlim([timestart timeend]);
 xlabel('Time (min)');
 
ax3=subplot(1,3,3);
% ax4=subplot(10,1,8:10);
 ft = edgeSensorVals(:,:,3);
% %rac1 = edgeSensorVals(:,:,1)./edgeSensorVals(:,:,2);
 imagesc(ft);
 set(gca, 'XTick',xlocs2, 'Xticklabel',xlabels2);
 set(gca, 'YTick',[], 'Yticklabel',[]);
% colormap(ax4,parula);
 climits = [prctile(vect(ft),5) prctile(vect(ft),95)];
caxis(ax3,climits);
 title('PLS3/F-tractin');
 xlabel('Time (min)');
 colorbar('EastOutside');
  
  xlim([timestart timeend]);


 
%%
%Visualize the masks/ratios
rangeX =[];

    rangeY = [];
    [imHight imWidth]=size(labelMask{1});
    TimeIndexes = idx_first:2:150;

    
    imRangeX = [];
    imRangeY = [];
    cellMontage = [];
    RGBFinalC1 = ones(imHight,imWidth);
     RGBFinalC2= ones(imHight,imWidth);
    RGBFinalC3 = ones(imHight,imWidth);
    for T = 1:length(TimeIndexes)
        time = TimeIndexes(T);

        thismask = (labelMask{time});
        %axis image
        wholeMask = thismask;
        wholeMask =imfill(wholeMask);
        perim = bwperim(wholeMask);
        perim = imdilate(perim, strel('disk',2));

    
    cmap = jet(length(TimeIndexes));
    cmap = [0,0,0; cmap];
    [imHight imWidth]=size(wholeMask);
    rgbIm=ones([imHight,imWidth,3]);

    perimVals = (perim)>0;
    RGBFinalC1(perimVals==1) = cmap(T,1);
    RGBFinalC2(perimVals==1) = cmap(T,2);
    RGBFinalC3(perimVals==1) = cmap(T,3);

    hold on;
    end
    RGBim(:,:,1) =  RGBFinalC1;
     RGBim(:,:,2) =  RGBFinalC2;
      RGBim(:,:,3) =  RGBFinalC3;
      figure(3);
      imshow(RGBim)
    %imagesc(cellMontage)
    axis image
    
 