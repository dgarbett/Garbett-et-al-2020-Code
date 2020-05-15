% Image processing for ladder patterns to track and idenitfy which pattern
% the cell is on
%
% Anjali Bisaria 2018
%


%% Working folder
clear
clc
close all
warning off

%what direction do file separations use 
filesep = '/'; %Anjali I think we don't need to define this, you can just use 'filesep' and it returns system default... right?

%Files
root = '/Volumes/AB_Data_3/DG_tplastin/20190920_BJ5 CRISPR analysis/WT'; %copy directory info into here
rawdir=[root,filesep];
datadir=[root,filesep,'data'];

prefix = 'Capture'; %prefix for loading files based on condition or prefix from 3i
output = 'Tracking_V1';%
moviebin = 2; %binning (1 or 2)

visualOutput = 1; %output some intermediate movies?
RGBwrite = 1; %needed for tracking code, keep at 1

runTracking = 1; %run cell tracking
segmentationCh = 0; %what are the channels used for segmentation; as in C0 or C1 etc
trackCh = 1; % what channel has the track information; as in C0 or C1 etc
timeStep = .5; %minutes between time steps;
numLanes = 1;




%make the appropriate directies
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
if moviebin == 2

    umPerPx = .325; %set to spatial resolution of camera/obj
    nucr = 20;
    debrisarea = 500;
    boulderarea = 200000;
    memdebrisarea=500;
    coff = 100;
    trackWidth = 20/umPerPx;
elseif moviebin == 1
    umPerPx = .1625;
    nucr = 80;
    debrisarea = 4000;
    boulderarea = 8000;
    coff = 200;
    memdebrisarea=5000;
    boulderarea=40000;
    trackWidth = 20/umPerPx;
end


close all



%% Execute Code
%

%read in the segmentation Channel
for ii = 1:length(segmentationCh)
    files = strcat('*',prefix, '*C',num2str(segmentationCh(1)),'*.tif');
    IMChStack{ii} = dir(files);
end


numFiles = numel(IMChStack{1,1});
%for each file make a tiff stack of the ratio;
for I = 1:numFiles
    
    
    
    IMChStackI = IMChStack{1,1}(I);
    IMChName{I} = IMChStackI.name;
    IMInfo{I} = imfinfo(IMChStackI.name);
    %what is the prefix for the filename from before channel info, assumes a multi-tiff stack with C0-CN
    filename = regexprep(IMChName{I},['_C',num2str(segmentationCh(1)),'.tif'],'');
    stripefile = regexprep(IMChName{I},['_C', num2str(segmentationCh(1)),'.tif'],['_C',num2str(trackCh(1)),'.tif']);
    
    stripeInfo = imfinfo(stripefile);
    %what frames to analyze, change if there are known problems with acquisition
    startFrame = 1;%assume start frame is 1 unless other specified
    endFrame = numel(IMInfo{1,I}); % image until end of stack
    
    frames = startFrame:endFrame;
    badframes = NaN(endFrame,1);
    
    %tracking parameters
    jitters = zeros(endFrame,2);
    jitterwindow=3*nucr;
    blobthreshold = -0.02;  %default 10xbin1=-0.02 20xbin2=-0.03;
    memblocksize = 500;
    maxmemnum = 100; % how many cells can we add max
    maxjumpmem = nucr; %how far can a cell move between frames to be tracked together

    
    
    
    if runTracking == 1
        %set up trace data files
        memtracedata = NaN(maxmemnum,endFrame,5);
        memtracking = NaN(maxmemnum,5);
        
        %set up initial tracking parameters form first 3 images        
        thisStack = IMChStack{1,1}(I);
        thisInfo = IMInfo{I};
        im_mem = double(imread(IMChName{I}, 1, 'Info', thisInfo));
        
        
        %tracking initialization
        [firstgoodindex,badframes,height,width] = timelapsesetup_PRIMO(im_mem,badframes);
        folder = [rawdir,filesep];
        
        %go through all frames
        for i = firstgoodindex:endFrame
            f = frames(i);
            
            thisInfo = IMInfo{I};
            %read in the segmentation channel
            rawMem = double(imread(IMChName{I}, f, 'Info', thisInfo));
            
            %segment
            mem_mask = segment_AdaptThresh(rawMem,debrisarea,51);
            %mem_mask=bwareaopen(mem_mask,debrisarea); %remove small objects
            mem_mask = bwperim(mem_mask); %just take the outline
            mem_mask = imfill(mem_mask,'holes'); %fill the outline 
            mem_mask = imerode(mem_mask,strel('disk',2)); %erode a bit
            mem_mask = imopen(mem_mask,strel('disk',1)); %get rid of small strings etc
            
            mem_mask  = imclearborder(mem_mask); %get rid of cells that touch the border
            mem_mask= bwareaopen(mem_mask,debrisarea); %get rid of small objects
            %mem_mask = imopen(mem_mask,strel('disk',3));
            
            [mem_label,nummem] =  bwlabel(mem_mask); %number cell objects
            
            stripeIM = double(imread(stripefile, f, 'Info',stripeInfo));
            stripe_mask = segment_AdaptThresh(stripeIM,5,51); %segment stripes
            

            stripeAngle = struct2cell(regionprops(stripe_mask,'Orientation')); %get the average "angle" of the stripes
            stripeAngle = (cell2mat(stripeAngle));
            stripeAngle(stripeAngle <0) = stripeAngle(stripeAngle <0) + 180;
            avstripeAngle = median(stripeAngle); %average across all segmented stripes
            stripe_fill = imclose(stripe_mask, strel('line', floor(trackWidth/2),avstripeAngle)); % expand the stripes given that angle
            stripe_fill = bwareaopen(stripe_fill,40); % get rid of small fragments
            stripe_label = bwlabel(stripe_fill); 

            mem_info = struct2cell(regionprops(mem_mask,rawMem,'Area','Centroid','MeanIntensity')'); %get info about each segmented cell
            mem_area = squeeze(cell2mat(mem_info(1,1,:)));
            mem_center = squeeze(cell2mat(mem_info(2,1,:)))';
            mem_density = squeeze(cell2mat(mem_info(3,1,:)));
            mem_mass = mem_density.*mem_area;
            
            %put information into memcurdata
            if numel(mem_area) == 0
                continue
            elseif numel(mem_area)== 1
                memcurdata = [mem_center(1),mem_center(2),mem_area,mem_mass'];
            else
                memcurdata = [mem_center(:,1),mem_center(:,2),mem_area,mem_mass];
            end
            
            
            if i>1
                lastgoodframe = find(badframes==0,1,'last'); %infd the last "good frame"
                H2BvsNLS = 1; %1:H2B 2:NLS ... not important for us
                
                %track based on last good frame
                [memtracedata,memcurdata,memtracking,mem_label] = adaptivetrack_membrane_PRIMO(lastgoodframe,f,memtracedata,memcurdata,memtracking,rawMem,mem_label,nucr,jitters(f,:),maxjumpmem,H2BvsNLS);
                badframes(f)=0;
            else
                badframes(f)=0;
            end
            
            %write the masks
            extractmemmask = bwmorph(mem_label,'remove');
            %mem_label_ring = mem_label - imerode(mem_label,strel('disk',1,0));
            
            %write output files of various masks used later calculations 
            if f == 1
                imwrite(uint16(mem_label),[RGBdir,filename,'_memedge.tif']); %outline of membranes labeled to match memtracedata
                imwrite(uint16(stripe_label),[RGBdir,filename,'_stripeedge.tif']); %imputed stripe outlines
                imwrite(uint16(stripe_mask),[RGBdir,filename,'_striperaw.tif']); %raw stripe outlines before the dilation
            else
                imwrite(uint16(mem_label),[RGBdir,filename,'_memedge.tif'],'WriteMode','append');
                imwrite(uint16(stripe_label),[RGBdir,filename,'_stripeedge.tif'],'WriteMode','append');
                imwrite(uint16(stripe_mask),[RGBdir,filename,'_striperaw.tif'],'WriteMode','append');
            end
            
            
            %cell and membrane tracking
            memid = find(~isnan(memcurdata(:,1)));
            numtrackmem = numel(memid);
            
            %get centroid and area information for membranes
            memcurdata = memcurdata(memid,:);
            mem_center = memcurdata(:,[1 2]); %center coordinates
            mem_area = memcurdata(:,3); %areas in px
            mem_mass = memcurdata(:,4); %total fluoresence 
            mem_info = regionprops(mem_label,'PixelIdxList'); 
            %mem_raw1 = regionprops(mem_label,imMem_raw,'PixelValues');
            mem_raw2 = regionprops(mem_label,rawMem,'PixelValues'); %get pixel values within cell masks
            mem_nuc_num =[];
            
            memIDvec = {mem_info.PixelIdxList};
            memidvec = ~cellfun(@isempty,memIDvec );
            
            mem_info = mem_info(memidvec);
            %mem_raw1 = mem_raw1(memidvec);
            mem_raw2 = mem_raw2(memidvec);
            %need to make sure only iterate through mem_info if it isn't empty
            if size(mem_center,1) == 1 && ~isempty(mem_raw2)
                mem_nuc_num(1,1) = mode(mem_raw2.PixelValues(mem_raw2.PixelValues>0)); %what number has been assigned to the membrane
            elseif ~isempty(mem_raw2)
                for m = 1:size(mem_center,1)
                    mem_nuc_num(m,1) = mode(mem_raw2(m).PixelValues(mem_raw2(m).PixelValues>0));
                end
            end;
            
            %add back to memtrace data array
            memtracedata(memid,f,:) = [mem_center,mem_area,mem_mass,mem_nuc_num];
            if maxmemnum - max(memid) < memblocksize
                memtempdata = NaN(memblocksize,endFrame,5);
                memtemptrack = NaN(memblocksize,5);
                memtracedata = [memtracedata;memtempdata];
                memtracking = [memtracking;memtemptrack];
                maxmemnum = maxmemnum+memblocksize;
            end
        end;
        
        %link up all the membranes 
        memtracedata(:,:,size(memtracedata,3)+1) = repmat([1:size(memtracedata,1)]',1,size(memtracedata,2));
        [memtracedata,~,~] = postprocessing_nolinking(memtracedata,memid,jitters,badframes,memtracking,maxmemnum,nucr);
        meminfo = imfinfo([RGBdir,filename,'_memedge.tif']);

        save([outputdir,'memtracedata_',filename,'.mat'],'memtracedata');
        
    end %runTracking ==1;
    % start ratio analysis
    close all;
end;
%%

%%%%%%%
%Section 2: individual cell analyses 

%read in all memtracedata files
files = dir('memtracedata*.mat');

for f = 1:numel(files)
    load(files(f).name);
    filename = regexprep(files(f).name,'memtracedata_','');
    filename = regexprep( filename ,'.mat','');
    cellNumber  = nanmean(memtracedata(:,:,6),2);
    %run through tracg ked cells
    for m = 1:length(cellNumber)
        cellNo = cellNumber(m);
        saveFile = 1; %switch in case a cell sucks and we don't want to output the data
        meminfo = imfinfo([RGBdir,filename,'_memedge.tif']);% load tracking tiff
        stripeinfo = imfinfo([RGBdir,filename,'_stripeedge.tif']); % load stripe edge tiff
        stripeoutlineinfo = imfinfo([RGBdir,filename,'_striperaw.tif']); % load  raw stripe tracking tiff
        cellData=[];
        
        
        cellData = squeeze(memtracedata(m,:,:));
        
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

         %calculate Mean Squared Displaceent from the first position
        MSD = pdist2([cellData(idx_first,1:2)],[cellData(idx_first:idx_last,1:2)])*umPerPx;
        
        %start looking at contact with stripes based on the # of tracked frames
        allStripes = [];
        stripe = {};
        numStripes  = nan(size(cellData,1),1);
        adhesionArea  = nan(size(cellData,1),1);
        for frames=1:length(frames_tracked) %iterate through tracked frames;
            frame = frames_tracked(frames);
            memEdgeImage=imread([RGBdir,filename,'_memedge.tif'],frame, 'Info', meminfo);
            memEdgeImageCell=uint16(memEdgeImage== cellNo); %only take segmented area of cell
            
            stripeImage = imread([RGBdir,filename,'_stripeedge.tif'],frame, 'Info', stripeinfo);
            uniqueStripes = stripeImage.*(memEdgeImageCell); %which stripes have overlap with membrane
            whichStripes = unique(uniqueStripes);
            whichStripes = whichStripes(2:end);
            stripe{frame} = whichStripes; %which stripes were occupied at this time
            allStripes = unique([allStripes ; whichStripes]);  %which stipes have ever been occupied
            numStripes(frame) = length(whichStripes); %how many stripes were occupied
            
            actualStripe = imread([RGBdir,filename,'_striperaw.tif'],frame, 'Info', stripeoutlineinfo);
            overlap = uint16(actualStripe>0) .* memEdgeImageCell;
            adhesionArea(frame) = sum(overlap(:)); %how much overlap between stripes and cell = adhesion area
        end
        
        %code to automatically detect which stripe the cells were located on
        stripeImage = imread([RGBdir,filename,'_stripeedge.tif'],idx_last, 'Info', stripeinfo);
        StripeInfo = regionprops(stripeImage, 'Centroid');
        stripeCenters = cell2mat(struct2cell(StripeInfo)');
        occupiedCenters = stripeCenters([whichStripes],:); 
        if size(occupiedCenters,1) > 1
            trackGap = mean(diff(occupiedCenters(:,1))*.325); %what is the averge distance between the occupied stripes
        else
            trackGap = 0;
        end;
        
        %assign the stripe
        if numLanes > 1
            clusterLanes = kmeans(stripeCenters(:,2),numLanes); %cluster stripes based on the vertical axis
            tbl = tabulate(clusterLanes);
            rowOrder = sortrows(tbl,2,'descend'); %the closest together gaps will have te highest number of stripes
            rowOrder = [rowOrder [1:numLanes]'];
            rowOrder = sortrows(rowOrder,1,'ascend');
            for l = 1:length(clusterLanes)
                clust = clusterLanes(l);
                LaneAssigned(l) = rowOrder(clust,4);
            end
            
            cellLane = mode(LaneAssigned([whichStripes]));
        else
            cellLane = 1; %if there are not different lane widths
        end
        
        
        %save ak te cell data
        cellData = [cellData numStripes adhesionArea];
        save([outputdir,filename,'_cell_analysis_',num2str(cellNo),'.mat'],'cellData','cellLane','trackGap','stripe','allStripes');
        
    end
end


%%
%combine all data

parentdir = '/Volumes/AB_Data_3/DG_tplastin/20190920_BJ5 CRISPR analysis';

minutesTracked = 60; %how long to check for tracking to keep a cell.  keep this the same to compare between replciates/conditions
subdirlist = {'WT' 'KO'}; %subdirs to look for cell data
trackingName = 'Tracking_V1'; 
cd(parentdir);
folderlist = dir(parentdir);
SummaryTable = {};
counter = 0;
timeStep = .5; %how long between imaging
umperpx = .325;
for s = 1:numel(subdirlist)
   % try
        disp(['Condition is ' subdirlist{s}]);
        cd([parentdir  '/' subdirlist{s} '/' trackingName]);
        allCellAnalysis = dir('*cell_analysis*.mat'); %read all cell data
        for a = 1:numel(allCellAnalysis)
            load(allCellAnalysis(a).name)
                counter = counter+1;
                longestrun = regionprops(~isnan(cellData(:,1)),'Area','PixelIdxList'); %look through longest tracked 
                if isempty(longestrun)
                    continue;
                end
                [sortedAreas, sortIndexes] = sort([longestrun.Area], 'Descend');
                idx_first = longestrun(sortIndexes(1)).PixelIdxList(1);
                idx_last = longestrun(sortIndexes(1)).PixelIdxList(end);
                
                %WT or mutant? (based on subdirlist)
                SummaryTable{counter,1} = s;
                %which exp did it come from 
                SummaryTable{counter,2} = allCellAnalysis(a).name;
                
                  
                [~, endIndex] = regexp(allCellAnalysis(a).name,'cell_analysis_');
                [startIndex,~] = regexp(allCellAnalysis(a).name,'.mat');
                %get the cell number
                cellNumber = extractBetween(allCellAnalysis(a).name,endIndex+1, startIndex-1);
                
                %cellnumber
                SummaryTable{counter,3} = str2num(cellNumber{1});
                
                %which Track is it on
                SummaryTable{counter,4} = cellLane;
                
    
                %when is it on a stripe?
                longeststripe = regionprops(~isnan(cellData(:,3)) & cellData(:,7)>0 ,'Area','PixelIdxList');
                if isempty(longeststripe)
                    continue;
                end
                
                [sortedAreas, sortIndexes] = sort([longeststripe.Area], 'Descend');
                stripe_first = longeststripe(sortIndexes(1)).PixelIdxList(1);
                stripe_last = longeststripe(sortIndexes(1)).PixelIdxList(end);

              
                
                %starting cell area
                SummaryTable{counter,5} = cellData(stripe_first,3)*umperpx^2;
                
                timeTracked = stripe_last-stripe_first;
                SummaryTable{counter,6} =timeTracked;
                numhrs = (timeTracked*timeStep);
               
                firstTimeVar = 6;
                if (timeTracked*timeStep) >= minutesTracked
                    %dArea after 1 hour
                    SummaryTable{counter,firstTimeVar+1} = (cellData(stripe_first+(minutesTracked/timeStep),3) - cellData(stripe_first,3))*umperpx^2;
                    %increase in number of bars after 1 hour
                    SummaryTable{counter,firstTimeVar+2}=  (cellData(stripe_first+(minutesTracked/timeStep),7) - cellData(stripe_first,7));
                else
                    SummaryTable{counter,firstTimeVar+1} = nan;
                    SummaryTable{counter,firstTimeVar+2} = nan;
                end

        end
    end
%end % through conditions
%end % throught folders

%Summary table variables: [1:subdir(#) 2:filename 3:cell(#) 4:cell lane 5:starting area in um2 6:time tracked (frames) 7:area change after 1 hour 8:change in bars after 1 hour

%%
for i = 1:size(SummaryTable,1)
    for j = 1 :size(SummaryTable,2)
        if isempty(SummaryTable{i,j})
            SummaryTable{i,j} = NaN; %fil empty values with nan
        end
        
        
    end
    % if sum(isnan(SummaryTable{i,5}))
end
%%
TableAll = (cell2table(SummaryTable));
varnames = {'Cell_Line','Cell_Folder','Cell_Number','Cell_Lane','Start_Area','tracked_time','hrdArea','hrdStripes'};
TableAll.Properties.VariableNames = varnames;
%%
%what are the cells you want to analyze
singleCell = cell2mat(SummaryTable(:,5)) < 5000; %initial tracked area needs to be smaller than this
%singlePost =   cell2mat(SummaryTable(:,8)) < 3; %can set to make sure only sat down on a single post 
Tracked = cell2mat(SummaryTable(:,6)) > minutesTracked/timeStep;
celltype = cell2mat(SummaryTable(:,1)) ==  1 %set depending on what subdir you want to analuze
ForAnalysis = TableAll(singleCell&Tracked&celltype,:);

allVals = NaN(1,size(ForAnalysis,1),size(ForAnalysis,2) -4);

hrarea = cell2mat(ForAnalysis.hrdArea);
starea = cell2mat(ForAnalysis.Start_Area);
stripenumber = cell2mat(ForAnalysis.Initial_Stripes);
hrarea(~isnan(hrarea)&starea<5000& stripenumber < 3);

%%
