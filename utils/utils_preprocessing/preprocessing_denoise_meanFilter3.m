%Use 3D mean filter to denoise input image sequence
%%ELiiiiiii, 20210916
function inputStack = preprocessing_denoise_meanFilter3(inputStack,bacthSizeArray)
%% inputs:
%     inputStack: 3D stack
%     bacthSizeArray: median filter size. If not set, this function
%                automatically sets it to [3,3,1].
%% check inputs
if nargin<2
    bacthSizeArray = [3,3,1];
end
if ~isvector(bacthSizeArray)
    error('Input bacthSizeArray must be a vector!');
end
if length(bacthSizeArray) ~= 3
    if length(bacthSizeArray) > 3
        warning('Too many elements in bacthSizeArray... Use the first 3');
        bacthSizeArray = bacthSizeArray(1:3);
    elseif length(bacthSizeArray) == 2
        warning('Not enough elements in bacthSizeArray...');
        bacthSizeArray = [bacthSizeArray(1),bacthSizeArray(1),bacthSizeArray(2)];
    elseif length(bacthSizeArray) == 1
        warning('Not enough elements in bacthSizeArray...');
        bacthSizeArray = [bacthSizeArray(1),bacthSizeArray(1),1];
    end    
end


%% 3d mean filter
inputStack = imboxfilt3(inputStack,bacthSizeArray);
end