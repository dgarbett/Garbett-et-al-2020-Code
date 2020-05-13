function cytolbl = segmentImageUsingThreshAndSeparate(cyto,varargin)
% General two step segmentation of a cytoplasmic (usually) marker image
% First step uses os thresholding, and the second tries to identify markers
% that are inward and uses some methd to separate the cells based on the
% inner marker
%
% segments an image using the strategy of dual segmentation
% first doing GMM type thresholding that gives a generous threshold
% that using this as a mask, use otsu to thrsold inner cell
% than using simple watershed I could identify the seperation between cells
% and apply that to the original gm threshold

arg.filename='';
arg.thresh='gm' ;%'gm';%'gm';
arg.separate='InnerTheshAndWatershed';
arg.minnucarea=50;
arg.mincellarea=500;
arg.suppressoversegmentation=1;
arg.redo=false;
arg.gaussfilter=fspecial('gaussian',[25 25],12);
arg.hmax=0.01;
arg=parseVarargin(varargin,arg);

minCellArea=arg.mincellarea;
minNucArea=arg.minnucarea;
suppressOverSegmentationParameter=arg.suppressoversegmentation;

%% if file exist - reread from drive
if ~isempty(arg.filename) && exist(arg.filename,'file') && ~arg.redo
    fprintf('reading file %s from disk\n',arg.filename)
    cytolbl = imread(arg.filename);
    return
end

%% background correct image 
cytoflat = backgroundSubtraction(cyto);

%% threshold 
switch arg.thresh
    case 'otsu'
        cytobw = optThreshold(cytoflat,'method','otsu');    
    case 'gm'
        cytobw = optThreshold(cytoflat,'method','gm','subsample',ceil(numel(cytoflat)*0.1));
    case 'logotsu'
        cytobw = optThreshold(cytoflat,'method','otsu','transform','log');
    otherwise
        error('unsupported thresholding method')
end
cytobw = bwareaopen(cytobw,minCellArea);

%% get inner segmentation 
switch arg.separate
    case 'InnerTheshAndWatershed'
        %% identif seeds by inner thresholding 
        appnucbw = optThreshold(cytoflat,'method','otsu','msk',cytobw);
        appnucbw = bwareaopen(appnucbw,minNucArea);
        dst = bwdist(~appnucbw);
        dst = imhmax(dst,suppressOverSegmentationParameter);
        wtrshd = watershed(-dst);
        cytolbl = bwlabel(cytobw & wtrshd>0);
        cytolbl = bwlabel(bwareaopen(cytolbl>0,minCellArea));
        
    case 'filterAndRegionamMaxAndNN'
        %% identify max and use watershed
        cyto2=imfilter(cytoflat,arg.gaussfilter);
        bwmx=imregionalmax(cyto2) & cytobw;
        dst = bwdist(bwmx);
        wtrshd = watershed(dst);
        cytolbl = bwlabel(cytobw & wtrshd>0);
        cytolbl = bwlabel(bwareaopen(cytolbl>0,minCellArea));
    case 'filterAndRegionamMaxAndWatershed'
        %% identify max and use watershed
        cyto2=imfilter(cytoflat,arg.gaussfilter);
        bwmx=imextendedmax(cyto2,arg.hmax);
        dst = bwdist(~bwmx);
        dst = imhmax(dst,suppressOverSegmentationParameter);
        wtrshd = watershed(-dst);
        cytolbl = bwlabel(cytobw & wtrshd>0);
%         cytolbl = bwlabel(bwareaopen(cytolbl>0,minCellArea));
    otherwise
        cytolbl =imfill(cytobw,'holes'); %AH: fills holes
end
%% if asked for write files
if ~isempty(arg.filename)
    imwrite(uint16(cytolbl),arg.filename);
    if ~isempty(regexp(arg.filename,'lbl','once'))
        imwrite(lbl2rgb(cytolbl),regexprep(arg.filename,'lbl','rgb'));
    end
end