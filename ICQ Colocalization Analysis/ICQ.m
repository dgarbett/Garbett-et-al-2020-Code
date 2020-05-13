% ICQ.m calculates the intensity correlation quotient.
% Inputs:
% Image_r/g-The raw images
% Red/Green_mu-The mean intensities (should be global mean)
% Sig_Map_r/g-The binary signal maps for the two channels
% Nuc_Map-the logical map identifying either the regions of the image to
% analyze

% Outputs:
% icqCoeff - the ICQ coefficient

% Updated 10 January 2019. Created 8/28/2018 By Dannielle G. McCarthy - Moerner Lab
function [icqCoeff] = ICQ(Image_r, Image_g, Red_mu, Green_mu, Sig_Map_r, Sig_Map_g, Nuc_Map)

%looks at all pixels in the image
if ~exist('Nuc_Map', 'var') && ~exist('Sig_Map_r', 'var') && ~exist('Sig_Map_g', 'var')
    
    %calculates the factors in the pdm
    red = Image_r - Red_mu;
    green = Image_g - Green_mu;
    
    %looks at all pixels in a specified region
elseif exist('Nuc_Map', 'var') && ~exist('Sig_Map_r', 'var') && ~exist('Sig_Map_g', 'var')
    
    %selects pixels within the specified region
    ind_nuc = find(Nuc_Map>0);
    
    %calculates the factors in the pdm
    red = Image_r(ind_nuc) - Red_mu;
    green = Image_g(ind_nuc) - Green_mu;
    
    %looks only at pixels with either red or green signal
else
    
    %selects pixels where there is either red or green signal
    ind_nuc_r = find(Sig_Map_r>0);
    ind_nuc_g = find(Sig_Map_g>0);
    ind_nuc = union(ind_nuc_r, ind_nuc_g);
    
    %calculates the factors in the pdm
    red = Image_r(ind_nuc) - Red_mu;
    green = Image_g(ind_nuc) - Green_mu;
    
end

%calculates the ICQ coeff
pdm = red.*green;

%counts the instances where pdm == 0 and both pixels are equal to the
%mean
mean_r = find(red(:)==0);
mean_g = find(green(:)==0);
eq_mean = intersect(mean_r, mean_g);

%counts the number of instances where the two channels varied in the
%same way
N_pos = sum(pdm(:) > 0);

%adds in instances where both channels were equal to the mean
N_same = N_pos + numel(eq_mean);

N_notsame = sum(pdm(:) <= 0) - numel(eq_mean);
%calculates the total number of pixels counted
N_tot = N_same + N_notsame;

icqCoeff = N_same/N_tot - 0.5;

end