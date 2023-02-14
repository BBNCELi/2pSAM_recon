% A simple function to expand a matrix
% Can also be implemented with padarray
%
% ELi, 20230213, add comments
function outputMatrix=expandAndCentralizeAMatrix(originalMatrix,outputMatrix_siz,fillinNumber)
%% size
[siz1,siz2]=size(originalMatrix);
if siz1>outputMatrix_siz(1) || siz2>outputMatrix_siz(2)
    error('No need to expand the matrix!');
end
%% centralize the input
outputMatrix = fillinNumber .* ones(outputMatrix_siz(1),outputMatrix_siz(2));
outputMatrix(round(outputMatrix_siz(1)/2)-floor(siz1/2):round(outputMatrix_siz(1)/2)+round(siz1/2)-1,round(outputMatrix_siz(2)/2)-floor(siz2/2):round(outputMatrix_siz(2)/2)+round(siz2/2)-1) = originalMatrix;
end