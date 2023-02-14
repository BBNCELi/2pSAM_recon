% A simple function to perform fftn for stacks
%
% NOTES
%     Using fftshift(fftn(stack)) directly assumes that the coordinate
%     locates at stack(0,0,0), which is not usually the case. ifftshift before
%     fftn equals to shifting the coordinate to stack(middle,middle,middle),
%     no matter what the stack size is.
%
% ELi, 20230213, add comments
function output = fftForStack(input)

output = fftshift(fftn(ifftshift(input)));