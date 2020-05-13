% sing_dist.m calculates the colocalization metric for a randomization of inputted images
% Inputs:
% meth - specifies the colocalization metric to be used with 0 for Pearson, 1 for ICQ, and 2 for Manders
% x/y_psf - denotes the x/y size of a PSF for scrambling in the images
% Image_r/g-The raw images
% Red/Green_mu-The mean intensities (should be global mean)
% Nuc_Map-the logical map identifying the regions to scramble
% im_scram - indicates which images you should scramble (0 for red, 1 for green, 2 for both, and 3 for neither)

% Outputs:
% coloc - outputs the parameter for measuring colocalization of the scrambled images
% scram_r/g - the scrambled image

%relies on code Manders.m, Pearson.m, ICQ.m

% Updated 13 February 2019. Created 8/28/2018 By Dannielle G. McCarthy - Moerner Lab

function [coloc, scram_r, scram_g] = sing_dist(meth, im_scram, x_psf, y_psf, Image_r, Image_g, Red_mu, Green_mu, Nuc_Map, Sig_Map_r, Sig_Map_g, red_thresh, green_thresh)

%assigns the signal thresholds to a deflaut value if needed
if ~exist('red_thresh', 'var') && ~exist('green_thresh', 'var')
    red_thresh = 0;
    green_thresh = 0;
end

%checks to see if the image(s) should scrambled
if im_scram ~= 3
    %assigns the output of the scrambled images
    scram_r = Image_r;
    scram_g = Image_g;
    
    if ismember(2, meth) && ~exist('Sig_Map_r', 'var') && ~exists('Sig_Map_g', 'var')
        
        error('must specify signal maps for Manders analysis');
        
    else
        
        if ~exist('Nuc_Map', 'var') && ~exist('Sig_Map_r', 'var') && ~exist('Sig_Map_g', 'var')
            
            %checks to see if a simple scramble is needed
            if x_psf == 1 && y_psf ==1
                
                %generates the linear indices for the image(assuming they
                %are the same size)
                [x_leng, y_leng] = size(scram_r)
                lin_ind = linspace(1, x_leng.*y_leng, x_leng.*y_leng);
                
                %generates a permutation of the linear indices and
                %scrambles the images
                if im_scram == 0
                    
                    rand_r = randperm(x_leng.*y_leng);
                    scram_r(lin_ind) = scram_r(rand_r);
                    
                elseif im_scram == 1
                    
                    rand_g = randperm(x_leng.*y_leng);
                    scram_g(lin_ind) = scram_g(rand_g);
                    
                else
                    rand_r = randperm(x_leng.*y_leng);
                    rand_g = randperm(x_leng.*y_leng);
                    scram_r(lin_ind) = scram_r(rand_r);
                    scram_g(lin_ind) = scram_g(rand_g);
                end
                
            else
                
                %computes the size the image assuming the images from each channel have
                %the same size
                [y_leng, x_leng] = size(scram_r);
                
                %computes the number of full PSF lengths along x = columns (and the
                %remainder)
                x_P = fix(x_leng./x_psf);
                x_rem = rem(x_leng, x_psf);
                
                %computes the number of full PSF lengths along y = rows (and the
                %remainder)
                y_P = fix(y_leng./y_psf);
                y_rem = rem(y_leng, y_psf);
                
                %generates the dimension vectors for the cell matrix conversion
                x_dim = ones(1, x_P).*x_psf;
                y_dim = ones(1, y_P).*y_psf;
                
                %appends on the last <psf-sized value if needed
                %MAY ONLY WORK IF THE PSF-SIZE IS DIVISBLE BY THE DIMENSIONS OF THE
                %IMAGE!!!
                if x_rem ~= 0
                    x_dim = [x_dim x_rem];
                end
                
                if y_rem ~= 0
                    y_dim = [y_dim y_rem];
                end
                
                %converts the image to a cell matrix in prep for scrambling
                cell_r = mat2cell(Image_r, y_dim, x_dim);
                cell_g = mat2cell(Image_g, y_dim, x_dim);
                
                %calculates the size of the cell matrix
                [cr, cc] = size(cell_r);
                
                %generates linear indices values for the images
                lin_ind = linspace(1, cr.*cc, cr.*cc);
                
                %calculates a different random permutation for the red and green
                %cell matrix indices and scrambles the images
                if im_scram == 0
                    
                    rand_r = randperm(cr*cc);
                    cell_r(lin_ind) = cell_r(rand_r);
                    
                elseif im_scram == 1
                    
                    rand_g = randperm(cr*cc);
                    cell_g(lin_ind) = cell_g(rand_g);
                    
                elseif im_scram ==2
                    
                    rand_r = randperm(cr*cc);
                    rand_g = randperm(cr*cc);
                    
                    %scrambles the cell matrix and converts back to matrix form
                    cell_r(lin_ind) = cell_r(rand_r);
                    cell_g(lin_ind) = cell_g(rand_g);
                    
                else
                    error('wrong syntax for indicating the image(s) to scramble')
                end
                
                %converts back to a matrix
                scram_r = cell2mat(cell_r);
                scram_g = cell2mat(cell_g);
                
            end
        else
            
            %checks to see if just a simple scramble is needed
            if x_psf == 1 && y_psf == 1
                
                %selects pixels within the nucleus
                ind_nuc = find(Nuc_Map>0);
                
                %computes the number of pixels in the nucleus
                num_pixels = numel(ind_nuc);
                
                %generates a permutation of the linear indices and
                %scrambles the images
                if im_scram == 0
                    
                    rand_r = randperm(num_pixels);
                    shuffle_r = ind_nuc(rand_r);
                    
                    %shuffles the pixels in the images
                    scram_r(ind_nuc) = scram_r(shuffle_r);
                    
                elseif im_scram == 1
                    
                    rand_g = randperm(num_pixels);
                    shuffle_g = ind_nuc(rand_g);
                    
                    %shuffles the psfs in the images
                    scram_g(ind_nuc) = scram_g(shuffle_g);
                    
                else
                    
                    rand_r = randperm(num_pixels);
                    rand_g = randperm(num_pixels);
                    
                    shuffle_r = ind_nuc(rand_r);
                    shuffle_g = ind_nuc(rand_g);
                    
                    %shuffles the psfs in the images
                    scram_r(ind_nuc) = scram_r(shuffle_r);
                    scram_g(ind_nuc) = scram_g(shuffle_g);
                    
                end
                
            else
                %computes the size of the image assuming the images from each channel have
                %the same size
                [y_leng, x_leng] = size(scram_r);
                
                %scrambles all regions at once
                
                %computes the number of full PSF lengths along x = columns (and the
                %remainder)
                x_P = fix(x_leng./x_psf);
                x_rem = rem(x_leng, x_psf);
                
                %computes the number of full PSF lengths along y = rows (and the
                %remainder)
                y_P = fix(y_leng./y_psf);
                y_rem = rem(y_leng, y_psf);
                
                %generates the dimension vectors for the cell matrix conversion
                x_dim = ones(x_P, 1).*x_psf;
                y_dim = ones(y_P, 1).*y_psf;
                
                %converts the images and binary maps to a cell matrix in prep for
                %scrambling excluding the edges if needed
                cell_r = mat2cell(Image_r(1:(y_leng - y_rem), 1:(x_leng - x_rem)), y_dim, x_dim);
                cell_g = mat2cell(Image_g(1:(y_leng - y_rem), 1:(x_leng - x_rem)), y_dim, x_dim);
                cell_nuc = mat2cell(sum(Nuc_Map(1:(y_leng - y_rem), 1:(x_leng - x_rem)),3), y_dim, x_dim);
                
                %calculates the size of the cell matrix
                [cr, cc] = size(cell_r);
                
                %reshapes the array to speed up the loop
                nuc_find = reshape(cell_nuc, [cr*cc, 1]);
                
                %matrix of the linear indices of psf completely inside of the
                %nucleus
                ind_nuc = zeros(cr*cc, 1);
                
                %loop to find the indices of the cells that are completely
                %within the nucleus
                for m = 1:cr*cc
                    
                    if all(all(nuc_find{m, 1}))
                        
                        ind_nuc(m) = m;
                        
                    end
                    
                end
                
                %deletes all of the zero entries from ind_nuc
                ind_nuc(ind_nuc==0)= [];
                
                %computes the number of pixels in the nucleus
                num_pixels = size(ind_nuc, 1);
                
                %calculates a different random permutation for the red and green
                %cell matrix indices and scrambles the images
                if im_scram == 0
                    
                    rand_r = randperm(num_pixels);
                    shuffle_r = ind_nuc(rand_r);
                    
                    %shuffles the psfs in the images
                    cell_r(ind_nuc) = cell_r(shuffle_r);
                    
                    
                elseif im_scram == 1
                    
                    rand_g = randperm(num_pixels);
                    shuffle_g = ind_nuc(rand_g);
                    
                    %shuffles the psfs in the images
                    cell_g(ind_nuc) = cell_g(shuffle_g);
                    
                elseif im_scram ==2
                    
                    rand_r = randperm(num_pixels);
                    rand_g = randperm(num_pixels);
                    
                    shuffle_r = ind_nuc(rand_r);
                    shuffle_g = ind_nuc(rand_g);
                    
                    %shuffles the psfs in the images
                    cell_r(ind_nuc) = cell_r(shuffle_r);
                    cell_g(ind_nuc) = cell_g(shuffle_g);
                    
                else
                    error('wrong syntax for indicating the image(s) to scramble')
                end
                
                %converts back to a matrix
                scram_r = cell2mat(cell_r);
                scram_g = cell2mat(cell_g);
                
                %appends the unscrambled parts of the original images to
                %the scrambled images if needed
                if x_rem ~= 0 && y_rem == 0
                    
                    scram_r = [scram_r Image_r(:, (x_leng - x_rem + 1): x_leng)];
                    scram_g = [scram_g Image_g(:, (x_leng - x_rem + 1): x_leng)];
                    
                elseif y_rem ~= 0 && x_rem == 0
                    
                    scram_r = [scram_r; Image_r((y_leng - y_rem + 1):y_leng, :)];
                    scram_g = [scram_g; Image_g((y_leng - y_rem + 1):y_leng, :)];
                    
                elseif x_rem ~= 0 && y_rem ~= 0
                    
                    %add the additional columns
                    scram_r = [scram_r Image_r(1:(y_leng - y_rem), (x_leng - x_rem + 1): x_leng)];
                    scram_g = [scram_g Image_g(1:(y_leng - y_rem), (x_leng - x_rem + 1): x_leng)];
                    
                    %adds the additional rows
                    scram_r = [scram_r; Image_r((y_leng - y_rem + 1):y_leng, :)];
                    scram_g = [scram_g; Image_g((y_leng - y_rem + 1):y_leng, :)];
                    
                end
            end
        end
    end
else
    %assigns the output of the scrambled images
    scram_r = Image_r;
    scram_g = Image_g;
end

%calculates the size of meth to determine if multiple analysis methods
%should be used
num_anal = numel(meth);
curr_anal = 1;

for j = 1:num_anal
    %performs the appropriate correlation/colocalization analysis on the
    %scrambled images
    
    this_meth = meth(j);
    
    if ~exist('Nuc_Map', 'var') && ~exist('Sig_Map_r', 'var') && ~exist('Sig_Map_g', 'var')
        
        if this_meth == 0
            
            coloc(j) = Pearson(scram_r, scram_g, Red_mu, Green_mu);
            
        elseif this_meth == 1
            
            coloc(j) = ICQ(scram_r, scram_g, Red_mu, Green_mu);
            
        elseif this_meth == 2
            
            error('signal map necessary for Manders analysis')
            
        else
            error('analysis method not specified by user inputs')
        end
        
    else
        
        %does not calculate new signal maps if the images were not scrambled
        if im_scram == 3
            
            if this_meth == 0
                
                coloc(curr_anal) = Pearson(scram_r, scram_g, Red_mu, Green_mu, Nuc_Map, Sig_Map_r, Sig_Map_g);
                curr_anal = curr_anal + 1;
                
            elseif this_meth == 1
                
                coloc(curr_anal) = ICQ(scram_r, scram_g, Red_mu, Green_mu, Sig_Map_r, Sig_Map_g, Nuc_Map);
                curr_anal = curr_anal + 1;
                
            elseif this_meth == 2
                
                [Mand_r, Mand_g, Frac_r, Frac_g, Frac_sig] = Manders(Sig_Map_r, Sig_Map_g, scram_r, scram_g, Nuc_Map);
                
                coloc(curr_anal:curr_anal + 4) = [Mand_r, Mand_g, Frac_r, Frac_g, Frac_sig];
                curr_anal = curr_anal + 5;
                
            else
                error('analysis method not specified by user inputs')
            end
        else
            %calculates new signal maps for the scrambled image(s)
            scram_Sig_Map_r = Sig_finder(scram_r, red_thresh, Nuc_Map);
            scram_Sig_Map_g = Sig_finder(scram_g, green_thresh, Nuc_Map);
            
            %calculates the colocalization coefficient using the new signal
            %maps
            if this_meth == 0
                
                coloc(curr_anal) = Pearson(scram_r, scram_g, Red_mu, Green_mu, Nuc_Map, scram_Sig_Map_r, scram_Sig_Map_g);
                curr_anal = curr_anal + 1;
                
            elseif this_meth == 1
                
                coloc(curr_anal) = ICQ(scram_r, scram_g, Red_mu, Green_mu, scram_Sig_Map_r, scram_Sig_Map_g, Nuc_Map);
                curr_anal = curr_anal + 1;
                
            elseif this_meth == 2
                
                [Mand_r, Mand_g, Frac_r, Frac_g, Frac_sig] = Manders(scram_Sig_Map_r, scram_Sig_Map_g, scram_r, scram_g, Nuc_Map);
                
                coloc(curr_anal:curr_anal + 4) = [Mand_r, Mand_g, Frac_r, Frac_g, Frac_sig];
                curr_anal = curr_anal + 5;
                
            else
                error('analysis method not specified by user inputs')
            end
            
        end
    end
end

end