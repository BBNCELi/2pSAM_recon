% Check inputs, load existed PSFs or generate new ones
% ELi, 20220530
function [psf_all,PSFParameters] = PSFsGenerator(...
    PSFParameters,                                             ...
    systemParameters_path,                                     ...
    saveLoc,                                                   ...
    NARatio,                                                   ...
    pinholeR_DL,                                               ...
    angDistx_normed, angDisty_normed,                          ...
    pupilPhase                                                 ...
    )
% PSFParameters: a struct containing parameters used when generating PSFs
%     .imagingMode: estimated NARatio and pinholeR_DL of 2pSAM system
%                   'mid' / 'min' . If not given, turn to function inputs "NARatio" "pinholeR_DL"
%     .EIParameters: parameters used when imaging samples (EI: elementary image)
%         .zoomFactor: zoom factor used. Determines FOV and pixel size in both PSFs and EIs. Default: 6
%         .catchSize: size of EIs. Default: 512
%         .upSample: upsampling rate in deconvolution. Determines FOV and pixel size in PSFs and EIs. Default: 1
%     .xysize: xysize of PSFs. should be LARGE ENOUGH to contain PSFs with large defocus and tilted angles
%         Default: PSFParameters.EIParameters.catchSize * PSFParameters.EIParameters.upSample
%     .angDistParameters: the distribution of angles for scanning
%         .scanningMode: angle-scanning mode: 0 / 1 / 2 / 3
%             %%switch
%             1: square scanning in back pupil plane
%         .SamplingNumber_X: sampling number x
%         .SamplingNumber_Y: sampling number y
%             2 (default): circle scanning in back pupil plane
%         .lowNA_rouNumber: how many circles in total. Default 2
%         .lowNA_thetaNumberMin: angle number in the inner circle. Default 4
%         .lowNA_thetaNumberMax: angle number in the outer circle. Default 8
%             3: sunflower scanning in back pupil plane
%         .NANums: how many angles in total
%             OTHERWISE: turn to function inputs "angDistx_normed" "angDisty_normed" to determine the anlge distribution
%             then set .scanningMode = 0
%             %%%%%%%
%         .Interval_X_Upper: max X. Default: 8.5
%         .Interval_Y_Upper: max Y. Default: 8.5
%         .bias_X: bias x. Default: 0
%         .bias_Y: bias y. Default: 0
%     .defocusUpper: max z (um) for PSF. Default: 50
%     .defocusLower: min z (um) for PSF. Default: -50
%     .defocusInterval: z interval (um) for PSF. Default: 1. defocusLower:defocusInterval:defocusUpper
%     .pupilPhaseSet: phase in pupil when generating PSFs
%         %%switch
%         'NO': zero-phase (default)
%         'RA': random phase. Followed by a number to set magnification
%         'SP': spherical aberration. Followed by a number to set magnification
%         'SE': estimated system phase stored in [systemParameters_path,'//pupilPhase//sePhase.mat']. Followed by a number to set magnification
%         ESPECIALLY: if function input "pupilPhase" is given, use it and set .pupilPhaseSet to timenow.
%         %%%%%%%
%     .energySet: axial energy normalization of PSFs. 0/1/2/3/4. Default: 2 (normalized for each angle)
% systemParameters_path: system parameters location. Containing:
%     //mid
%     //min
%     //pupilPhase
%     //angleDistInPupil
%     overall.mat
% saveLoc: PSFs save location

%% check inputs
% overall
load([systemParameters_path,'//overall.mat'],'numAper','lambda_nm','rindexObj','rindexSp','FOV_zoom8_um');
    PSFParameters.numAper=numAper;
    PSFParameters.rindexObj=rindexObj;
    PSFParameters.rindexSp=rindexSp;
    PSFParameters.lambda=lambda_nm;
    PSFParameters.FOV_zoom8_um = FOV_zoom8_um;
% EIParameters
if ~isfield(PSFParameters,'EIParameters'); PSFParameters.EIParameters = struct(); end
if ~isfield(PSFParameters.EIParameters,'zoomFactor'); PSFParameters.EIParameters.zoomFactor = 6; end
if ~isfield(PSFParameters.EIParameters,'catchSize'); PSFParameters.EIParameters.catchSize = 128; end
if ~isfield(PSFParameters.EIParameters,'upSample'); PSFParameters.EIParameters.upSample = 1; end
% xysize
if ~isfield(PSFParameters,'xysize'); PSFParameters.xysize = PSFParameters.EIParameters.catchSize * PSFParameters.EIParameters.upSample; end
% imagingMode, pinhole size && NA ratio
if isfield(PSFParameters,'imagingMode') && (isequal(PSFParameters.imagingMode,'mid') || isequal(PSFParameters.imagingMode,'min'))
    if exist('pinholeR_DL','var') && ~isempty(pinholeR_DL)
        warning('Ignoring input: pinholeR_DL...');
    end
    if exist('NARatio','var') && ~isempty(NARatio)
        warning('Ignoring input: NARatio...');
    end
    pinholeR_DL = loadFirstVariable([systemParameters_path,'//',PSFParameters.imagingMode,'//pinholeR_DL.mat']);
    PSFParameters.pinholeR_DL = pinholeR_DL;
    NARatio = loadFirstVariable([systemParameters_path,'//',PSFParameters.imagingMode,'//NARatio.mat']);
    PSFParameters.NARatio = NARatio;
else
    if ~exist('pinholeR_DL','var') || isempty(pinholeR_DL)
        error('Unrecognized pinholeR_DL');
    end
    if ~exist('NARatio','var') || isempty(NARatio)
        error('Unrecognized NA ratio');
    end
    PSFParameters.pinholeR_DL = pinholeR_DL;
    PSFParameters.NARatio = NARatio;
    PSFParameters.imagingMode = ['pinhole',num2str(pinholeR_DL),'NA',num2str(NARatio)];
end
% angle number and their distribution
if isfield(PSFParameters,'angDistParameters') && ~isequal(PSFParameters.angDistParameters.scanningMode,0)
    if (exist('angDistx_normed','var') && ~isempty(angDistx_normed)) || (exist('angDisty_normed','var') && ~isempty(angDisty_normed))
        warning('Ignoring input: angDistx_normed and angDisty_normed...');
    end
    angDist_path = [systemParameters_path,'//angleDistInPupil'];
    if ~isfield(PSFParameters.angDistParameters,'Interval_X_Upper'); PSFParameters.angDistParameters.Interval_X_Upper = 8.5; end %xMax used when scanning
    if ~isfield(PSFParameters.angDistParameters,'Interval_Y_Upper'); PSFParameters.angDistParameters.Interval_Y_Upper = 8.5; end %yMax used when scanning
    if ~isfield(PSFParameters.angDistParameters,'bias_X'); PSFParameters.angDistParameters.bias_X = 0; end %xBias used when scanning
    if ~isfield(PSFParameters.angDistParameters,'bias_Y'); PSFParameters.angDistParameters.bias_Y = 0; end %yBias used when scanning
    [angDistx_normed, angDisty_normed] = angDistGenerator(PSFParameters,angDist_path);
    PSFParameters.angDistx_normed = angDistx_normed;
    PSFParameters.angDisty_normed = angDisty_normed;
else
    if (~exist('angDistx_normed','var') || isempty(angDistx_normed)) || (~exist('angDisty_normed','var') || isempty(angDisty_normed))
%         error('Unrecognized angle distribution');
        warning('Using default settings for angle-scanning...');
        PSFParameters.angDistParameters.scanningMode = 2; %0:customized; 1:square; 2:circle; 3:sunflower
        PSFParameters.angDistParameters.lowNA_rouNumber = 2;
        PSFParameters.angDistParameters.lowNA_thetaNumberMin = 4;
        PSFParameters.angDistParameters.lowNA_thetaNumberMax = 8;
        PSFParameters.angDistParameters.Interval_X_Upper = 8.5;
        PSFParameters.angDistParameters.Interval_Y_Upper = 8.5;
        PSFParameters.angDistParameters.bias_X = 0;
        PSFParameters.angDistParameters.bias_Y = 0;
        angDist_path = [systemParameters_path,'//angleDistInPupil'];
        [angDistx_normed, angDisty_normed] = angDistGenerator(PSFParameters,angDist_path);
    else
        PSFParameters.angDistParameters.scanningMode = 0;
    end
    PSFParameters.angDistx_normed = angDistx_normed;
    PSFParameters.angDisty_normed = angDisty_normed;
end
PSFParameters.angleNum = length(PSFParameters.angDistx_normed);
if ~isfield(PSFParameters,'defocusUpper'); PSFParameters.defocusUpper = 25; end
if ~isfield(PSFParameters,'defocusLower'); PSFParameters.defocusLower = -25; end
if ~isfield(PSFParameters,'defocusInterval'); PSFParameters.defocusInterval = 1; end
% phase in pupil, used to compensate aberrations
if ~(isfield(PSFParameters,'pupilPhaseSet')) || (exist('pupilPhase','var') && ~isempty(pupilPhase))
    if ~exist('pupilPhase','var') || isempty(pupilPhase)
        PSFParameters.pupilPhaseSet = 'NO';
    else
        PSFParameters.pupilPhaseSet = datestr(now, 'YYYYmmDD_HHMMSS');
        PSFParameters.duringGeneration.pupilPlaneMatrix = pupilPhase;
    end
end
% axial energy modulation
if ~(isfield(PSFParameters,'energySet')); PSFParameters.energySet = 0; end
% K
if ~(isfield(PSFParameters,'shiftKRange')); PSFParameters.shiftKRange = 0; end
if ~(isfield(PSFParameters,'shiftKBias')); PSFParameters.shiftKBias = 0; end
% systemParameters_path
PSFParameters.systemParameters_path = systemParameters_path;

%% save path
PSFParameters.PSFFolderName = PSFNameGenerator(PSFParameters);
savePath = [saveLoc,'//',PSFParameters.PSFFolderName];

%% load or generate
overwriteFlag = 0;
if exist(savePath,'dir')
    if overwriteFlag
        warning(['PSF in ',savePath,' will be overwrited...']);
        [psf_all,PSFParameters] = PSFsSimuGPU(PSFParameters,savePath);
    else
        disp('Found existed PSF. Loading now...');
        load([savePath,'/psf_all.mat'],'psf_all');
        load([savePath,'/PSFParameters.mat'],'PSFParameters');
    end
else
    mkdir(savePath);
    [psf_all,PSFParameters] = PSFsSimuGPU(PSFParameters,savePath);
end