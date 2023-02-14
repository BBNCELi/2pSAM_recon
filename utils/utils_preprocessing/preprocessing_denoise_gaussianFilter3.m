%Use 3D gaussian filter to denoise input image sequence
%%ELiiiiiii, 20211101
function inputStack = preprocessing_denoise_gaussianFilter3(inputStack,sigmaArray)
%% inputs:
%     inputStack: 3D stack
%     sigmaArray: sigma value of 3D gaussian filter. If not set, this function
%                automatically sets it to [3,3,1].

%% check inputs
if nargin<2
    sigmaArray = [3,3,1];
end
if ~isvector(sigmaArray)
    error('Input bacthSizeArray must be a vector!');
end
if length(sigmaArray) ~= 3
    if length(sigmaArray) > 3
        warning('Too many elements in bacthSizeArray... Use the first 3');
        sigmaArray = sigmaArray(1:3);
    elseif length(sigmaArray) == 2
        warning('Not enough elements in bacthSizeArray...');
        sigmaArray = [sigmaArray(1),sigmaArray(1),sigmaArray(2)];
    elseif length(sigmaArray) == 1
        warning('Not enough elements in bacthSizeArray...');
        sigmaArray = [sigmaArray(1),sigmaArray(1),1];
    end    
end

%% 3d gaussian filter
inputStack = imgaussfilt3(inputStack,sigmaArray,'padding','replicate');
end