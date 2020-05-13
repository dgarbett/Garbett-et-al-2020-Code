function bg_img=calculate_bg_img_rm_blobs(filepath,filenamecontains)
% Calculates an averaged background image based on multiple files containing
% the string "filenamecontains" located in the folder "path". Eliminates
% "dirt" on empty image frames. 
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
% Arnold Hayer, 5 June 2014

files=getFilenames(filepath);
files=files(boolRegExp(files,filenamecontains));

imstack_raw=[];
for i=1:length(files)
    tempim_raw=double(imread([filepath,filesep,files{i}]));
    tempim_filt1=imfilter(tempim_raw,fspecial('disk',5),'symmetric');
    imstack(:,:,i)=tempim_filt1;  
    disp(num2str(i));
end
imstack_mean=mean(imstack,3);
stack_NaN=[];

for i=1:length(files)
    tempframe=imstack(:,:,i);
    diff=(tempframe-imstack_mean);
    tempframe(diff>(0.5*std(tempframe(:))))=NaN;
    stack_NaN(:,:,i)=tempframe;
    disp(num2str(i));
end
bg_img=uint16(nanmean(stack_NaN,3));
bg_img=imfilter(bg_img,fspecial('disk',50),'symmetric');
% subplot(1,2,1),imagesc(imstack_mean);
% subplot(1,2,2),imagesc(bg_img);  
% imagesc(bg_img);
imwrite(bg_img,[filepath,filesep,'AVG_bg',filenamecontains,'.tif'],'tif','Compression','none');

