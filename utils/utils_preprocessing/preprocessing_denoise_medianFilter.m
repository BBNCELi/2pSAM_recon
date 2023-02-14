%Use median filter to denoise input image sequence
%%ELiiiiiii, 20210628
function inputStack = preprocessing_denoise_medianFilter(inputStack,bacthSize)
%% inputs:
%     inputStack: 3D stack
%     batchSize: median filter size. If not set, this function
%                automatically sets it to 3.
%% check inputs
if nargin<2
    bacthSize = 3;
end

%% median filter
for zCount = 1:size(inputStack,3)
    inputStack(:,:,zCount) = medfilt2(inputStack(:,:,zCount),[bacthSize,bacthSize]);
end

end