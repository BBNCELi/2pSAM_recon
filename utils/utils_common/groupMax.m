% This function group-max-projects an input stack according to given
% groupSize, cropFlag and projDim.
% if projDim is not given, this function projects inputStack along the 3rd(z) dimension.
% if cropFlag is not given, this funciton crops the inputStack.

%%ELiiiiiii, 20210819
function projResult = groupMax(inputStack, groupSize, cropFlag, projDim)
%% defaults
if nargin < 4
    projDim = 3;
end
if nargin < 3
    cropFlag = true;
end
if ndims(inputStack)>3
    error('The inputStack must be 3D');
end
if ~(projDim <= 3) || ~(projDim > 0) || ~(fix(projDim) == projDim)
    error('Incorrect project dimension');
end
[size1, size2, size3] = size(inputStack);

%% sizes
if cropFlag == true || rem(size3, groupSize) == 0
    projResultSize3 = floor(size3 / groupSize);
else
    if projDim ~= 3
        error('Can not perform a 1(x)/2(y)-dim projection while avoiding cropping!');
    end
    projResultSize3 = floor(size3 / groupSize) + 1;
end

if projDim == 3
    projResultSize1 = size1;
    projResultSize2 = size2;
elseif projDim == 1
    projResultSize1 = size2;
    projResultSize2 = groupSize;
elseif projDim == 2
    projResultSize1 = size1;
    projResultSize2 = groupSize;
end

projResult = zeros(projResultSize1, projResultSize2, projResultSize3, 'like', inputStack);

%% project group by group
for projCount = 1:projResultSize3
    projResult(:,:,projCount) = squeeze(max(inputStack(:,:,(projCount-1)*groupSize+1:min(projCount*groupSize,size3)),[],projDim));
end