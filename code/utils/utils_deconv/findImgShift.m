function [shifty,shiftx] = findImgShift(img_ref,img_move)

if size(img_ref)~=size(img_move)
    error('');
end

%%
% [yi,xi]=size(img_ref);
% 
% FI1 = fft2(img_ref);
% FI2 = fft2(img_move);
% FR = FI1.*conj(FI2);%calculating correlation 
% R = ifft2(FR);
% R = fftshift(R);
% 
% num = find(R==max(R(:)));
% 
% if length(num)~=1
%     num=num(round(length(num)/2));
% end
% 
% [i,j] = ind2sub(size(R), num);
% shifty = floor(yi/2)-i+1;
% shiftx = floor(xi/2)-j+1;
%%



%% it seems code below works better
c = gather(normxcorr2(img_move,img_ref));
[ypeak,xpeak] = find(c==max(c(:)));
if length(ypeak)~=1
    disp('Warning: inaccurate motion estimation!');
    ypeak=ypeak(round(length(ypeak)/2));
    xpeak=xpeak(round(length(xpeak)/2));
end
shifty = size(img_ref,1)-ypeak;
shiftx = size(img_ref,2)-xpeak;
%%