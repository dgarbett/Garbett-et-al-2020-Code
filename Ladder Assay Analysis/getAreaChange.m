parentdir = '/Volumes/AB_Data_3/DG_tplastin/1_9 area replot';
%filesdir = 'New Controls';
minutesTracked = 60;
subdirlist = {'WT'};
trackingName = 'Tracking_V1';
cd(parentdir);
folderlist = dir();
area = [];
counter = 0;
timeStep = .5;
umperpx = .325;
counter2 = 0;
for f= 4
    cd(parentdir);
    cd(folderlist(f).name);

        try

            allCellAnalysis = dir('*cell_analysis*.mat');
            for a = 1:numel(allCellAnalysis)
                load( allCellAnalysis(a).name)
                if isnan(cellLane)
                    continue
                else
                    counter2 = counter2 + 1;
                    
                    longestrun = regionprops(~isnan(cellData(:,1)),'Area','PixelIdxList');
                    if isempty(longestrun)
                        continue;
                    end
                    [sortedAreas, sortIndexes] = sort([longestrun.Area], 'Descend');
                    idx_first = longestrun(sortIndexes(1)).PixelIdxList(1);
                    idx_last = longestrun(sortIndexes(1)).PixelIdxList(end);
                    
                    
                    %when is it on a stripe?
                    longeststripe = regionprops(~isnan(cellData(:,3)) & cellData(:,7)>0 ,'Area','PixelIdxList');
                    
                    if isempty(longeststripe)
                        continue;
                    end
                    
                    %find the longest tracked run for the biggest object found for a cell (in case
                    %something small is left over after erosion etc)
                    [sortedAreas, sortIndexes] = sort([longeststripe.Area], 'Descend');
                    stripe_first = longeststripe(sortIndexes(1)).PixelIdxList(1);
                    stripe_last = longeststripe(sortIndexes(1)).PixelIdxList(end);
                    
                    
                    %starting cell area is single cell
                    if cellData(stripe_first,3)*umperpx^2 > 5000
                        continue;
                    end
                    %SummaryTable{counter,5} = cellData(stripe_first,3)*umperpx^2;
                    
                    
                    timeTracked = stripe_last-stripe_first;
                    
                    numhrs = (timeTracked*timeStep);
                    
                    if (timeTracked*timeStep) >= minutesTracked
                        counter = counter + 1;
                        numpts = minutesTracked/timeStep;
                        cellArea = cellData(:,3)*(umperpx^2);
                        area(counter,:) = cellArea(stripe_first:stripe_first+numpts);
                    else
                        continue;
                    end
                    
                    
                    
                end%on a known stripe
            end %through files
        end %try
   % end % through conditions
end % throught folders

%%
