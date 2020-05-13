% coloc.m processes all images specified by fig_titles and calculates the colocalization distribution from
% scrambling one image in each image pair.

% Inputs: (hard coded)
% sig_r/g-arrays of the desired signal thresholds
% user-defined inputs

% Outputs: (saved to a text file)
% coloc_icq - a 2D array with regions in columns and ICQ coefficient
%           at each signal threshold specified in sig_r/g
% dist_icq - a 2D array with regions in columns and ICQ coefficient
%            at each signal threshold specified in sig_r/g for the scrambling
%            specified by im_scram

%relies on code cell_mask.m, sing_dist.m

% Adapted 10 November 2019 By Dannielle G. McCarthy - Moerner Lab

function coloc()

fig_titles = {'Actin_T_plastin_lamellipodia', 'Actin_T_plastin_cytoplasm'};

%initializes .mat file to save the full ICQ data
save('Actin_T_plastin.mat', 'fig_titles');

sigrs = {[1], [1]};
siggs = {[1], [1]};

%user-input: which images should you scramble
im_scram = 1;%input('Which images should you scramble? (0 for red, 1 for green, 2 for both, and 3 for neither)');

%user-input: size of pixel blocks to scramble
x_psf = 1;%input('number of integer columns of the PSF');
y_psf = 1;%input('number of integer rows of the PSF');

%user-input: mean or median threshhold? (O for mean and 1 for median)
mean_med = 0;%input('use mean or median? (O for mean and 1 for median)');

% meth - specifies the colocalization metric to be used with 0 for Pearson, 1 for ICQ, and 2 for Manders
meth = 1;

[~,total_types] = size(fig_titles);

%user-input: subtract background using rolling ball technique?
sub_back = 0;%input('subtract background using rolling ball technique? (0 for no; 1 for yes)');

[~,total_types] = size(fig_titles);

for curr = 1:total_types
    %selects next image type
    fig_title = fig_titles{curr};
    %creates a text file to save data
    fileID = fopen(strcat(fig_title, '.txt'), 'w');
    %user-input: directory of images and file save format
    folder_path = 'C:\Users\dgmccart\Documents\MATLAB\Actin_collaboration';
    %image naming convention (in folder green images listed first, then red images)
    im_names = '\20190805_HUVEC_*STED*_z1_DTCorr.tiff';
    
    %creates a directory of images
    images = dir(strcat(folder_path, im_names));
    
    %writes general info to the text file
    fprintf(fileID, strcat('\r\n\r\n','data files\r\n'));
    
    %determines the number of images
    num = size(images);
    numIm = num(1);
    
    %total image sets
    numb = numIm./2; %2 for just green/red, 3 for 2 green, 1 red, 4 for STED only bkgd
    
    %counter for the number of images
    curr_im = 1;
    
    %loops over all the image sets
    for curr_set = 1:numb
        
        %loads first green image
        Im = images(curr_im);
        Image_g = double(imread(strcat(Im.folder, '\', Im.name)));
        
        %adds the image name to the text file
        fprintf(fileID, strcat(Im.name, '\r\n'));
        
        %subtract background using rolling ball technique
        if sub_back
            [Image_g, bkgd_g] = RollBall(Image_g, 30);
        end
        
        %Adds the green channel image to the green channel images output matrix
        Greens{curr_set, 1} = Image_g;
        
        %increments the image counter
        curr_im = curr_im + 1;
        
        %loads red image
        Im = images(curr_im);
        Image_r = double(imread(strcat(Im.folder, '\', Im.name)));
        
        %adds the image name to the text file
        fprintf(fileID, strcat(Im.name, '\r\n'));
        
        %subtract background using rolling ball technique
        if sub_back
            [Image_r, bkgd_r] = RollBall(Image_r, 30);
        end
        
        %Adds the red channel image to the red channel images output matrix
        Reds{curr_set, 1} = Image_r;
        
        %increments the image counter
        curr_im = curr_im + 1;
        
        %creates a mask of the are in the cell to analyze
        Cell_Map = cell_mask(Image_g);
        
        %Adds the Cell_Map image to the Cell_Map images output matrix
        Cell_Maps{curr_set, 1} = Cell_Map;
    end
    
    %save variables to a matlab file
    save(strcat(fig_title, '_images.mat'), 'Reds', 'Greens', 'Cell_Maps');
    
    %clear variables
    clear Reds Greens Cell_Maps
end

load(strcat(fig_title, '_images.mat'));
%sets colocalization outputs to zero to avoid errors
coloc_icq = [];
dist_icq = [];

%selects the appropriate signal thresholds for the two channels
sig_r = sigrs{curr};
sig_g = siggs{curr};

%selects next image type
fig_title = fig_titles{curr};
%creates a text file to save data
fileID = fopen(strcat(fig_title, '.txt'), 'w');

%writes general info to the text file
fprintf(fileID, strcat('\r\n\r\n', 'scrambled image(s)(0 for red, 1 for green, and 2 for both)  '));
fprintf(fileID,'%f',im_scram);
fprintf(fileID, strcat('\r\n\r\n', 'size of psf(x_psf, y_psf)'));
fprintf(fileID,'%f %f',[x_psf, y_psf]);

load(strcat(fig_title, '_images.mat'));

fprintf(fileID, strcat('\r\n\r\n', 'red channel signal threshold', '\r\n\r\n'));
fprintf(fileID,'%f\r\n', sig_r);
fprintf(fileID, strcat('\r\n\r\n', 'green channel signal threshold', '\r\n\r\n'));
fprintf(fileID,'%f\r\n', sig_g);

%calculates the size of the signal threshold input arrays
[num_sig, ~] = size(sig_r);

%calculates the number of images
[num_im, ~] = size(Cell_Maps);

for thresh = 1:num_sig
    
    %selects the signal thresholds for the 2 channels
    sig_thresh_r = sig_r(thresh);
    sig_thresh_g = sig_g(thresh);
    
    for ims = 1:num_im
        
        %selects the images and associated nuclear map
        curr_red = Reds{ims};
        curr_gr = Greens{ims};
        curr_nm = Cell_Maps{ims};
        
        %calculates signal maps considering all the nuclei at once
        
        Sig_Map_r = Sig_finder(curr_red, sig_thresh_r, curr_nm);
        Sig_Map_g = Sig_finder(curr_gr, sig_thresh_g, curr_nm);
        
        %adds the signal maps to the signal map 3D matrix
        sig_num = ims + (thresh - 1).*num_im;
        Sig_Maps_r{sig_num, 1} = Sig_Map_r;
        Sig_Maps_g{sig_num, 1} = Sig_Map_g;
        
        
        %finds the linear indices of pixels with signal in either the red or the
        %green channel
        ind_nuc_r = find(Sig_Map_r>0);
        ind_nuc_g = find(Sig_Map_g>0);
        ind_nuc = union(ind_nuc_r, ind_nuc_g);
        
        if mean_med == 0
            
            %calculates the mean signal in each channel
            Red_mu = round(mean2(curr_red(ind_nuc)));
            Green_mu = round(mean2(curr_gr(ind_nuc)));
            
        elseif mean_med == 1
            
            %calculates the median signal in each channel
            Red_mu = round(median(curr_red(ind_nuc)));
            Green_mu = round(median(curr_gr(ind_nuc)));
            
        else
            error('wrong syntax for specifying mean/median')
        end
        
        %adds mu or mean values to a matrix
        mu_red(ims, 1) = Red_mu;
        mu_green(ims, 1) = Green_mu;
        
        %calculates the colocalization coefficient on the original
        %images(im_scram = 3 --> do not scramble either images)
        [col_og, ~, ~] = sing_dist(meth, 3, x_psf, y_psf, curr_red, curr_gr, Red_mu, Green_mu, curr_nm, Sig_Map_r, Sig_Map_g, sig_thresh_r, sig_thresh_g);
        
        %calculates the colocalization coefficient on the appropriately scrambled images
        [col_rand, ~, ~] = sing_dist(meth, im_scram, x_psf, y_psf, curr_red, curr_gr, Red_mu, Green_mu, curr_nm, Sig_Map_r, Sig_Map_g, sig_thresh_r, sig_thresh_g);
        
        %adds value to the original colocalization distribution matrix
        curr_coloc_icq(ims, 1) = col_og(meth);
        
        %adds value to the randomly scrambled distribution matrix
        curr_dist_icq(ims, 1) = col_rand(meth);
        
        
        
    end
end

%save variables to a matlab file
save(strcat(fig_title, '_coloc.mat'), 'Sig_Maps_r', 'Sig_Maps_g');

%writes to the text file

%writes a header for the values
fprintf(fileID, strcat('\r\n\r\n', 'signal threshold in the red channel'));
fprintf(fileID,'%f', sig_thresh_r);
fprintf(fileID, strcat('\r\n\r\n', 'signal threshold in the green channel'));
fprintf(fileID,'%f', sig_thresh_g);

if mean_med == 0
    fprintf(fileID, strcat('\r\n\r\n', 'mean red signal', '\r\n\r\n'));
else
    fprintf(fileID, strcat('\r\n\r\n', 'median red signal', '\r\n\r\n'));
end
fprintf(fileID,'%f\r\n', mu_red);

if mean_med == 0
    fprintf(fileID, strcat('\r\n\r\n', 'mean green signal', '\r\n\r\n'));
else
    fprintf(fileID, strcat('\r\n\r\n', 'median green signal', '\r\n\r\n'));
end
fprintf(fileID,'%f\r\n', mu_green);

%adds value to the original colocalization distribution matrix
coloc_icq(thresh, :) = curr_coloc_icq;

%adds value to the randomly scrambled distribution matrix
dist_icq(thresh, :) = curr_dist_icq;

%write to the text file
fprintf(fileID, strcat('\r\n\r\n', 'ICQ colocalization coefficient', '\r\n\r\n'));
fprintf(fileID,'%f\r\n',coloc_icq);
fprintf(fileID, strcat('\r\n\r\n', 'ICQ colocalization coefficient after scrambling', '\r\n\r\n'));
fprintf(fileID,'%f\r\n',dist_icq);

save(strcat(fig_title, '_coloc.mat'),'coloc_icq', 'dist_icq');
clear Reds Greens Cell_Maps Sig_Maps_r Sig_Maps_g mu_red mu_green curr_coloc_icq curr_dist_icq
close all
fclose all;

end
