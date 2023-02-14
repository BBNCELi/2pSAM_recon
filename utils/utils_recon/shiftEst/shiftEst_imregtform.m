% Sub-pixel shift estimation using function "imregtform" provided by Matlab
% 
% NOTES
%     This method provides the highest sub-pixel accuracy at the cost of
%     speed.
% 
% Structured by Eli, 20221231
function [yDelta , xDelta] = shiftEst_imregtform(imgRef,imgMov)
%% check inputs
if ~ismatrix(imgRef)
    error('Shift estimation only supports images as inputs');
end
if ~isequal(size(imgRef), size(imgMov))
    error('The input images for shift estimation should have the same size');
end

%% make sure that the inputs for imregtform are not gpu variables
if isgpuarray(imgRef) || isgpuarray(imgMov)
    imgRef = gather(imgRef);
    imgMov = gather(imgMov);
end

%% imregtform
[optimizer,metric] = imregconfig('monomodal');
tform = imregtform(imgRef, imgMov, 'translation', optimizer, metric);
xDelta = tform.T(3,1);
yDelta = tform.T(3,2);

