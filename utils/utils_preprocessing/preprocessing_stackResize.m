% This function resamples the input 3-D images stack slice-by-slice.
% ELiiiiiii, 20210123.

function stack_after = preprocessing_stackResize(stack_before,sizeOfEachSlice)
% Inputs:
%     stack_before: a 3-D stack to be resampled.
%     sizeOfEachSlice: expected size of each slice in stack_after.
% Outputs:
%     stack_after: 3-D stack after resampling.

%% dimension judgement
dim_stack = ndims(stack_before);
switch(dim_stack)
    case 1
        warning('SUSPICIOUS: resampling an array using imresize...');
    case 2
        disp('Resampling an image...');
    case 3
        disp('Resampling a 3-D stack...');
    otherwise
        error('Resampling for 4-D or higher-dimension input is temporarily unsupported.');
end

%% preparations
stack_after_size = [sizeOfEachSlice,size(stack_before,3)];
stack_after = zeros(stack_after_size,'like',stack_before);

%% resize slice by slice
for sliceNum = 1:size(stack_before,3)
    stack_after(:,:,sliceNum) = imresize(stack_before(:,:,sliceNum),sizeOfEachSlice,'bilinear');
end