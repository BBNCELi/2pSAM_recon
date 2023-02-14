% Sub-pixel shift estimation using Foroosh2002 method
% 
% NOTES
%     Ref:
%     Foroosh H, Zerubia J B, Berthod M. Extension of phase correlation to subpixel registration [J]. 
%         Ieee Transactions on Image Processing, 2002, 11(3): 188-200.
%     Original code:
%     Lulu Ardiansyah (2022). ExtPhaseCorrelation 
%         (https://www.mathworks.com/matlabcentral/fileexchange/26417-extphasecorrelation), 
%         MATLAB Central File Exchange. Retrieved December 30, 2022.
% 
% Modified by Eli, 20221231
function [yDelta , xDelta] = shiftEst_Foroosh2002(imgRef, imgMov)
%% check inputs
if ~ismatrix(imgRef)
    error('Shift estimation only supports images as inputs');
end
if ~isequal(size(imgRef), size(imgMov))
    error('The input images for shift estimation should have the same size');
end

%% Phase correlation (Kuglin & Hines, 1975)
img1 = imgRef;
img2 = imgMov;
af = fftn(double(img1));
bf = fftn(double(img2));
cp = af .* conj(bf) ./ abs(af .* conj(bf));
cp(isnan(cp)) = 0;
cp = cp - mean(cp,'all'); % to remove the zero-frequency component after ifft
if ~any(cp,'all') % to avoid mistakes when isequal(imgRef,imgMov)
    yDelta = 0; xDelta = 0;
    return
end
icp = real(ifft2(cp));
mmax = max(max(icp));

%Find the main peak
[x,y,v] = find(mmax == icp);
%figure; imshow(icp,[],'notruesize');
if length(y) > 1
    disp('Multiple peaks are found... May lead to inaccurate shift estimation');
    y=y(round(length(y)/2));
    x=x(round(length(x)/2));
elseif isempty(y)
    disp('Found no peaks during shift estimation... Output zero shifts...');
    yDelta = 0; xDelta = 0;
    return
end

%% Extension to Subpixel Registration [Foroosh, Zerubia & Berthod, 2002]
[M, N] = size(img1);

%two side-peaks
xsp = x + 1;
xsn = x - 1;
ysp = y + 1;
ysn = y - 1;

%if the peaks outsize the image, then use xsn and/or ysn for the two
%side-peaks
if xsp > M
    xsp = M-1;
end
if ysp > N
    ysp = N-1;
end

%Calculate the translational shift
deltaX1 = ((icp(xsp,y) * xsp + icp(x,y) * x) / (icp(xsp,y) + icp(x,y)))-1;
deltaY1 = ((icp(x,ysp) * ysp + icp(x,y) * y) / (icp(x,ysp) + icp(x,y)))-1;

%I don't know why but after few test i find out that the result of deltaX
%and delta Y is inverted.. :( ??

%Validate if translation shift is negative
if deltaX1 < (N/2)
    deltaY = deltaX1;
else
    deltaY = deltaX1 - M;
end

if deltaY1 < (M/2)
    deltaX = deltaY1;
else
    deltaX = deltaY1 - N;
end


%% make sure that the outputs are not gpu variables for future convenience
yDelta = -deltaY;
xDelta = -deltaX;
if isgpuarray(yDelta)
    yDelta = gather(yDelta);
    xDelta = gather(xDelta);
end