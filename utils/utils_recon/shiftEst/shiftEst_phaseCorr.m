% Shift estimation using the classic phase correlation method
% 
% NOTES
%     Ref:
%     Kuglin C D, Hines D C. The phase correlation image alignment method [Z]. 
%         International Conference on Cybernetics and Society. 
%         San Francisco, CA, USA; IEEE, New York, NY, USA. 1975: 163-165
% 
% Implemented by Eli, 20221231
function [yDelta , xDelta] = shiftEst_phaseCorr(imgRef,imgMov)
%% check inputs
if ~ismatrix(imgRef)
    error('Shift estimation only supports images as inputs');
end
if ~isequal(size(imgRef), size(imgMov))
    error('The input images for shift estimation should have the same size');
end

%% the normalized cross power spectrum
imgRef_fft = fftshift(fftn(ifftshift(imgRef)));
imgMov_fft = fftshift(fftn(ifftshift(imgMov)));
pDif_fft = imgRef_fft .* conj(imgMov_fft) ./ abs(imgRef_fft .* conj(imgMov_fft));
pDif_fft(isnan(pDif_fft)) = 0;
pDif_fft = pDif_fft - mean(pDif_fft,'all'); % to remove the zero-frequency component after ifft
if ~any(pDif_fft,'all') % to avoid mistakes when isequal(imgRef,imgMov)
    yDelta = 0; xDelta = 0;
    return
end
pDif = real(fftshift(ifftn(ifftshift(pDif_fft))));
pDif_max = max(pDif, [], 'all');

%% find main peak
[yPeak, xPeak, ~] = find(pDif == pDif_max);
if length(yPeak) > 1
    disp('Multiple peaks are found... May lead to inaccurate shift estimation');
    yPeak=yPeak(round(length(yPeak)/2));
    xPeak=xPeak(round(length(xPeak)/2));
elseif isempty(yPeak)
    disp('Found no peaks during shift estimation... Output zero shifts...');
    yDelta = 0; xDelta = 0;
    return
end
[yDim, xDim] = size(imgRef);
yDelta = round((yDim+1)/2) - yPeak;
xDelta = round((xDim+1)/2) - xPeak;

%% make sure that the outputs are not gpu variables for future convenience
if isgpuarray(yDelta)
    yDelta = gather(yDelta);
    xDelta = gather(xDelta);
end