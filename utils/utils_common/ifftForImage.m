% A simple function to perform ifft2 for image
%
% NOTES
%     Using fftshift(ifft2(img)) directly assumes that the image coordinate
%     locates at img(0,0), which is not usually the case. ifftshift before
%     ifft2 equals to shifting the coordinate to img(middle,middle), no
%     matter what the image size is.
%
% ELi, 20230213, add comments
function output = ifftForImage(input)

if ~ismatrix(input)
    error('fftForImage supports only image input');
end

output = fftshift(ifft2(ifftshift(input)));