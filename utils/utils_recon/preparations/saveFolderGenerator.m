% Folder name generator for 2pSAM reconstruction
% 
% NOTES
%     According to name-generating rules set by ELi in 20210624
% 
%%ELi, 20210624
%%ELi, 20220531, add MC
%%ELi, 20230131, rearrange struct
function saveFolder = saveFolderGenerator(reconOpts,PSFParameters)
saveFolder = ['recon',...
    '_',...
    reconOpts.solver,...
    '_',...
    'zoom',num2str(PSFParameters.EIParameters.zoomFactor),...
    'ca',num2str(PSFParameters.EIParameters.catchSize),...
    'up',num2str(PSFParameters.EIParameters.upSample),...
    '_',...
    num2str(PSFParameters.defocusLower),'-',num2str(PSFParameters.defocusInterval),'-',num2str(PSFParameters.defocusUpper),...
    '_',...
    'pup',PSFParameters.pupilPhaseSet,...
    '_',...
    'en',num2str(PSFParameters.energySet),...
    '_',...
    'K',num2str(PSFParameters.shiftKRange),...
    'bias',num2str(PSFParameters.shiftKBias),...
    '_',...
    'DAO',num2str(reconOpts.DAO),...
    '_',...
    'EW',num2str(reconOpts.upWeight),...
    '_',...
    'mI',num2str(reconOpts.maxIter),...
    '_',...
    'mC',num2str(reconOpts.demotion),...
    '_',...
    'frames',num2str(reconOpts.startFrame),'-',num2str(reconOpts.step),'-',num2str(reconOpts.endFrame),...
    ];
end