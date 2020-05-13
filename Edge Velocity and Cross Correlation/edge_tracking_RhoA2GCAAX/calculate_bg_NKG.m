function bg_img=calculate_bg_NKG(filepath,filenamecontains)
% Calculates an averaged background image based on multiple files containing
% the string "filenamecontains" located in the folder "path". 
%
% Input 
%   path
%   filenamecontains
%
% Output 
%   averaged background image, writes the file into path
%
% Uses 
%   getFilenames.m
%   boolRegExp.m
%
% Arnold Hayer, 4 June 2015

subdir=getSubdirectories(filepath);
files=getFilenames(filepath);
files=files(boolRegExp(files,filenamecontains));

imstack_raw=[];
for i=1:length(subdir)
    files=getFilenames([filepath,filesep,subdir{i},filesep,'Pos0']);
    file=files(boolRegExp(files,filenamecontains));
    tempim_raw=double(imread([filepath,filesep,subdir{i},filesep,'Pos0',filesep,file{1}]));
%     tempim_filt1=imfilter(tempim_raw,fspecial('disk',5),'symmetric');
    imstack(:,:,i)=tempim_raw;  
    disp(num2str(i));
end

bg_img=uint16(mean(imstack,3));
imwrite(bg_img,[filepath,filesep,'AVG_bg',filenamecontains,'.tif'],'tif','Compression','none');

