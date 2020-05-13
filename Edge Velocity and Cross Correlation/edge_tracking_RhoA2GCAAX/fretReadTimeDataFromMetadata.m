function timeData=fretReadTimeDataFromMetadata(root,dataSubdir,numBefore)
%This function reads in the time data from the metadata files created by
%micromanager. It returns a cell array with one vector of numbers for each
%subdirectory (dataSubdir should be a cell array of subdirectory names) in
%the root folder specified by "root". A third input "numBefore" can be used
%to specify the number of frames before stimulus. Time zero will then be
%defined as the time half way in between frame numBefore and frame
%(numBefore + 1).

if nargin<3
    numBefore=0;
end

%% Read in time data
if ~isempty(dataSubdir)
    timeData=cell(size(dataSubdir));
    for i=1:length(dataSubdir)
        fid=fopen([root filesep dataSubdir{i} filesep 'metadata.txt'],'r');
        lines=readAllLines(fid);
        fclose(fid);
        frameKeyLines=find(boolRegExp(lines,'^"FrameKey'));
        timeLines=find(boolRegExp(lines,'^[\s]+"ElapsedTime-ms"'));
        timeLines=timeLines(timeLines>min(frameKeyLines));
        if length(frameKeyLines)~=length(timeLines)
            fprintf('Warning on subdir %s -- # of time lines does not match number of frame lines.\n',dataSubdir{i});
        end
        lines=lines(timeLines);
        tokens=getTokens(lines,'"ElapsedTime-ms"\: ([0-9]+),');
        dnums=str2double(tokens)/1000;
        if numBefore>0
            timeData{i}=(dnums-0.5*(dnums(numBefore)+dnums(numBefore+1)));
        else
            timeData{i}=dnums-dnums(1);
        end
    end
else
    folder=root;
    fid=fopen([folder filesep 'metadata.txt'],'r');
    lines=readAllLines(fid);
    fclose(fid);
    frameKeyLines=find(boolRegExp(lines,'^"FrameKey'));
    timeLines=find(boolRegExp(lines,'^[\s]+"ElapsedTime-ms"'));
    timeLines=timeLines(timeLines>min(frameKeyLines));
    if length(frameKeyLines)~=length(timeLines)
        fprintf('Warning on subdir %s -- # of time lines does not match number of frame lines.\n');
    end
    lines=lines(timeLines);
    tokens=getTokens(lines,'"ElapsedTime-ms"\: ([0-9]+),');
    dnums=str2double(tokens)/1000;
    if numBefore>0
        timeData=(dnums-0.5*(dnums(numBefore)+dnums(numBefore+1)));
    else
        timeData=dnums-dnums(1);
    end
end