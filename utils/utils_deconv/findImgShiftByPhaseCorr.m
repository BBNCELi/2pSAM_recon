%Implementation of:
%   Extension of Phase Correlation to Subpixel Registration
%   (Hassan Foroosh, Josiane B. Zerubia, Marc Berthod)
%Implemented by:
%   Lulu Ardiansyah (halluvme@gmail.com)
%
%TODO:
%   - Find out whether this implementation is correct :)
%   - Combine the result, overlay the images based on the result
%
%                                               eL-ardiansyah
%                                               January, 2010
%                                                       CMIIW
%============================================================
function [deltaY , deltaX] = findImgShiftByPhaseCorr(img1, img2)
%Description:
%   Find the translation shift between two image
%
%Parameter:
%   img1 = image 1
%   img2 = image 2
%   image 1 and image 2 , must in the same size

%Phase correlation (Kuglin & Hines, 1975)
%============================================================

%this estimation error is rather big compared to findImgShift WHEN the
%actually shifts are LARGE

%this method requires calculation precision during the fftn procedure
% i don't know why but if we don't transform img1 to a double, mmax will be
% NaN + NaN*i
img1_gather = double(gather(img1));
img2_gather = double(gather(img2));
af = fftshift(fft2(ifftshift(img1_gather)));
bf = fftshift(fft2(ifftshift(img2_gather)));
cp = af .* conj(bf) ./ abs(af .* conj(bf));
icp = real(fftshift(ifft2(ifftshift(cp))));
mmax = max(max(icp));

%Find the main peak
[y,x,v] = find(mmax == icp);
%figure; imshow(icp,[],'notruesize');

%Extension to Subpixel Registration [Foroosh, Zerubia & Berthod, 2002]
%============================================================
[M, N] = size(img1);

%two side-peaks
ysp = y + 1;
ysn = y - 1;
xsp = x + 1;
xsn = x - 1;

%if the peaks outsize the image, then use xsn and/or ysn for the two
%side-peaks
if ysp > M
    ysp = M-1;
end
if xsp > N
    xsp = N-1;
end

y_centralized = (M+1)/2 - y;
x_centralized = (N+1)/2 - x;
ysp_centralized = (M+1)/2 - ysp;
xsp_centralized = (N+1)/2 - xsp;

%Calculate the translational shift
deltaY = ((icp(ysp,x) * ysp_centralized + icp(y,x) * y_centralized) / (icp(ysp,x) + icp(y,x)));
deltaX = ((icp(y,xsp) * xsp_centralized + icp(y,x) * x_centralized) / (icp(y,xsp) + icp(y,x)));

%I don't know why but after few test i find out that the result of deltaX
%and delta Y is inverted.. :( ??

%Validate if translation shift is negative
% if deltaX1 < (N/2)
%     deltaY = -deltaX1;
% else
%     deltaY = M - deltaX1;
% end
% 
% if deltaY1 < (M/2)
%     deltaX = -deltaY1;
% else
%     deltaX = M - deltaY1;
% end