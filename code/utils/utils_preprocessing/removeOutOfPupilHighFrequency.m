function [conv_sum,imageParameters] = removeOutOfPupilHighFrequency(conv_sum,zoomFactor)
%ELiiiiiii, 20210611
conv_sum(conv_sum<0) = 0;
if size(conv_sum,1)~=size(conv_sum,2)
    error('XY size must be equal in this function');
end
imageParameters.xysize_zoomFactor = size(conv_sum,1);

%%
imageParameters.zoomFactor = zoomFactor;
imageParameters.dxy_zoomFactor = (86e3*8 / imageParameters.zoomFactor) / imageParameters.xysize_zoomFactor;%nm, dxy of conv_sum


imageParameters.numAper=1.05;
imageParameters.lambda=910;%nm
imageParameters.centrePixelIndex = round(imageParameters.xysize_zoomFactor/2);%PSFParameters.xysize/2+1;%(PSFParameters.xysize+1)/2;
imageParameters.dkxy = (2*pi)/(imageParameters.xysize_zoomFactor*imageParameters.dxy_zoomFactor);
imageParameters.kMax = (2*pi)/(2*(imageParameters.lambda/2/imageParameters.numAper));
imageParameters.kMaxPixels = imageParameters.kMax/imageParameters.dkxy;%radius

pupil=Circle([imageParameters.xysize_zoomFactor,imageParameters.xysize_zoomFactor],[imageParameters.centrePixelIndex,imageParameters.centrePixelIndex],[0,imageParameters.kMaxPixels]);

%%
for imageCount = 1:size(conv_sum,3)
    conv_sum(:,:,imageCount) = abs(ifftForImage(fftForImage(conv_sum(:,:,imageCount)).*pupil));
end

