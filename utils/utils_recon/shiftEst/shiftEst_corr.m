% Shift estimation using the simplest normalized-cross-correlation method
% 
% NOTES
%     This function serves both speed and pixel-level accuracy. It is also 
%     quite robust when large shifts are present
% 
% Implemented by Eli, 20221231
function [yDelta , xDelta] = shiftEst_corr(imgRef,imgMov)
%% check inputs
if ~ismatrix(imgRef)
    error('Shift estimation only supports images as inputs');
end
if ~isequal(size(imgRef), size(imgMov))
    error('The input images for shift estimation should have the same size');
end

%% correlation
c = normxcorr2(imgMov,imgRef);
[ypeak,xpeak] = find(c==max(c(:)));
if length(ypeak)~=1
    disp('Warning: inaccurate motion estimation!');
    ypeak=ypeak(round(length(ypeak)/2));
    xpeak=xpeak(round(length(xpeak)/2));
end
yDelta = size(imgRef,1)-ypeak;
xDelta = size(imgRef,2)-xpeak;

%% make sure that the outputs are not gpu variables for future convenience
if isgpuarray(yDelta)
    yDelta = gather(yDelta);
    xDelta = gather(xDelta);
end