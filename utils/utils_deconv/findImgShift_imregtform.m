function [deltaY , deltaX] = findImgShift_imregtform(img1, img2)

img1_gather = double(gather(img1));
img2_gather = double(gather(img2));

[optimizer,metric] = imregconfig('monomodal');
tform = imregtform(img1_gather, img2_gather, 'translation', optimizer, metric);
deltaX = tform.T(3,1);
deltaY = tform.T(3,2);