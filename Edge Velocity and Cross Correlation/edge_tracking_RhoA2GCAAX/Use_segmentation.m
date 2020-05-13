%% Segement_image.m
clear;clc;close all;
image=double(imread('test_image.tif'));
%image=imfilter(image_raw,fspecial('disk',1),'symmetric');

mask=segmentImageUsingThreshAndSeparate(image,'separate','');

mask_fine=imfilter(mask,fspecial('disk',2),'symmetric');
%mask_fine = refineMaskEdges(mask,image,3,5,0);

subplot(2,2,1);imshow(mat2gray(image));
subplot(2,2,2);imshow(mat2gray(mask));
subplot(2,2,3);imshow(mat2gray(mask_fine));

imRGB(:,:,1)=mat2gray(mask_fine);
imRGB(:,:,2)=mat2gray(image);
imRGB(:,:,3)=zeros(size(image));

subplot(2,2,4);imshow(imRGB);