% cell_mask.m - creates a binary mask from freehand drawn areas
% inputs:
% Image - matrix of image to make a mask for

% Adapted 12 Septmeber 2019 by DGM - Moerner Lab

function Cell_Map = cell_mask(Image)

%plots the image
figure
imagesc(Image);

% Ask user for the # of areas to mask
num_nuc = input('How many areas do you want to mask?');
close

Cell_Map = [];

for nuc = 1:num_nuc
    
    %logical to determine if a satisfactory cell map was made
    notdone = true;
    
    while notdone
        %draw a mask for the image
        figure
        imagesc(Image);
        
        h = imfreehand;
        partial_Cell_Map = h.createMask();
        close
        
        % Plot Image and binary map
        figure
        subplot(1,2,1)
        imagesc(Image)
        title(strcat('are you satisfied with the binary cell map for', num2str(nuc)));
        subplot(1,2,2)
        imagesc(partial_Cell_Map)
        pause
        close
        
        %user-input: continue code or remake map
        im_scram = input('Are you satisfied with the mask? (0 for yes and 1 for no)');
        
        if im_scram == 0
            notdone = false;
        end
        
    end
    
    Cell_Map(:, :, nuc) = partial_Cell_Map;
    
end

Cell_Map = sum(Cell_Map, 3);

if num_nuc ~=1
    % Plot Image and binary map
    figure
    subplot(1,2,1)
    imagesc(Image)
    subplot(1,2,2)
    imagesc(Cell_Map)
    pause
end
end
