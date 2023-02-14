% A simple function to perform ifftn for stacks
%
% NOTES
%     Using fftshift(ifftn(stack)) directly assumes that the coordinate
%     locates at stack(0,0,0), which is not usually the case. ifftshift before
%     ifftn equals to shifting the coordinate to stack(middle,middle,middle),
%     no matter what the stack size is.
%
% ELi, 20230213, add comments
function output = ifftForStack(input)

output = fftshift(ifftn(ifftshift(input)));