% Get updating sequence when doing multi-angle deconvolution
%ELi, 20220531, add central angle index in output
function [outputSequence,centreAngleIndex] = getIterSequence(PSFParameters)
%% preparations
angleNum = PSFParameters.angleNum;
angleScanningMode = PSFParameters.angDistGenerationMode;
shiftx = PSFParameters.angDistx_normed;
shifty = PSFParameters.angDisty_normed;
shiftxy = sqrt(shiftx.^2 + shifty.^2);
%% switch
switch angleScanningMode
    case 0 % unknown input angle distribution
        [~, outputSequence] = sort(shiftxy,'descend');
        outputSequence(shiftxy > 1) = [];
        centreAngleIndex = outputSequence(end);
    case 1 % square scanning
        [~, outputSequence] = sort(shiftxy,'descend');
        outputSequence(shiftxy > 1) = [];
        centreAngleIndex = outputSequence(end);
    case 2 % circle scanning
        outputSequence = [2:angleNum,1];
        outputSequence(shiftxy > 1) = [];
        centreAngleIndex = 1;
    case 3 % sunflower scanning
        [~, outputSequence] = sort(shiftxy,'descend');
        outputSequence(shiftxy > 1) = [];
        centreAngleIndex = outputSequence(end);
end
outputSequence = outputSequence(:)';
end
