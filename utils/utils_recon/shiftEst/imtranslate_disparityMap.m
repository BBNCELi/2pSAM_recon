% Translate image according to given disparity map
%
% NOTES
% 
% ELi, 20230208
function img = imtranslate_disparityMap(img, dispMap)
% INPUT
%     img       - input image, 2d matrix
%     dispMap   - disparity map, [img_r,img_c,2]

%% mesh grid
[img_r,img_c] = size(img);
[X,Y] = meshgrid(1:img_r,1:img_c);
Xq = X - dispMap(:,:,2);
Yq = Y - dispMap(:,:,1);

%% translate image
img = interp2(X,Y,img,Xq,Yq,'bilinear',0);