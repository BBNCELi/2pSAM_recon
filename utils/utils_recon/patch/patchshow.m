% A simple function to imshow image patches
% ELi, 20230104
function patchshow(I,contrast)
[r,c] = size(I);
if nargin < 2; contrast = [minInCell(I), maxInCell(I)]; end
hold on;
for rcount = 1:r
    for ccount = 1:c
        subplot(r,c,(rcount-1)*c+ccount);
        imshow(I{rcount,ccount},contrast);
    end
end
hold off;