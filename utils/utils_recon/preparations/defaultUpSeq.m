% Get updating sequence for volume reconstruction in 2pSAM
% 
% NOTES
%     This sequence has no impact to the final results but affect the 
%     convergence speed in a great sense. Quite similar to OS (ordered-
%     subset) reconstruction in X-ray CT.
% 
% REFERENCE
%     Beister M, Kolditz D, Kalender W A. Iterative reconstruction methods
%     in X-ray CT [J]. Physica Medica-European Journal of Medical Physics,
%     2012, 28(2): 94-108.
% 
% ELi, 20220531, add central angle index in output
% ELi, 20230209, add comments
function [upSeq,CAInd] = defaultUpSeq(PSFParameters)
% INPUT
%     PSFParameters - the parameters of PSFs used for imaging

%% extract for better use
angleNum = PSFParameters.angleNum;
angleScanningMode = PSFParameters.angDistGenerationMode;
shiftx = PSFParameters.angDistx_normed;
shifty = PSFParameters.angDisty_normed;
shiftxy = sqrt(shiftx.^2 + shifty.^2);

switch angleScanningMode
    %% unknown angle-scanning method
    case 0
        [~, upSeq] = sort(shiftxy,'descend');
        upSeq(shiftxy > 1) = [];
        CAInd = upSeq(end);
    %% square scanning
    case 1
        [~, upSeq] = sort(shiftxy,'descend');
        upSeq(shiftxy > 1) = [];
        CAInd = upSeq(end);
    %% circle scanning (default)
    case 2
        upSeq = [2:angleNum,1];
        upSeq(shiftxy > 1) = [];
        CAInd = 1;
    %% sunflower scanning
    case 3
        [~, upSeq] = sort(shiftxy,'descend');
        upSeq(shiftxy > 1) = [];
        CAInd = upSeq(end);
end
upSeq = upSeq(:)';

%% projections in the same cell unit will be update simultaneously
upSeq = num2cell(upSeq);
end
