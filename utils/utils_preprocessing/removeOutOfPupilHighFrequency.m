% Remove the frequency components outside pupil circle
%
% ELi, 20230214, add comments
function [proj_all,imageParameters] = removeOutOfPupilHighFrequency(proj_all,zoomFactor)
%% size
proj_all(proj_all<0) = 0;
if size(proj_all,1)~=size(proj_all,2)
    error('XY size must be equal in this function');
end
imageParameters.xysize_zoomFactor = size(proj_all,1);

%% defaul imaging parameters
imageParameters.zoomFactor = zoomFactor;
imageParameters.dxy_zoomFactor = (86e3*8 / imageParameters.zoomFactor) / imageParameters.xysize_zoomFactor;%nm, dxy of conv_sum
imageParameters.numAper=1.05;
imageParameters.lambda=920;%nm
imageParameters.centrePixelIndex = round(imageParameters.xysize_zoomFactor/2);%PSFParameters.xysize/2+1;%(PSFParameters.xysize+1)/2;
imageParameters.dkxy = (2*pi)/(imageParameters.xysize_zoomFactor*imageParameters.dxy_zoomFactor);
imageParameters.kMax = (2*pi)/(2*(imageParameters.lambda/2/imageParameters.numAper));
imageParameters.kMaxPixels = imageParameters.kMax/imageParameters.dkxy;%radius
pupil=Circle([imageParameters.xysize_zoomFactor,imageParameters.xysize_zoomFactor],[imageParameters.centrePixelIndex,imageParameters.centrePixelIndex],[0,imageParameters.kMaxPixels]);

%% filtering
for imageCount = 1:size(proj_all,3)
    proj_all(:,:,imageCount) = abs(ifftForImage(fftForImage(proj_all(:,:,imageCount)).*pupil));
end