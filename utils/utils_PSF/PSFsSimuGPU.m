% Prepare parameters, generate and save new PSFs
%%ELiiiiiii, 20220530
function [psf_all,PSFParameters] = PSFsSimuGPU(PSFParameters,savePath)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% preparations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basics
PSFParameters.EIParameters.zoomFactor;
PSFParameters.EIParameters.xysize = PSFParameters.EIParameters.catchSize*PSFParameters.EIParameters.upSample;
PSFParameters.dxy = (PSFParameters.FOV_zoom8_um*1e3*8 / PSFParameters.EIParameters.zoomFactor) / PSFParameters.EIParameters.xysize;%nm, EQUAL TO pixel size of EI
PSFParameters.xysize;
PSFParameters.centrePixelIndex = round(PSFParameters.xysize/2);
PSFParameters.duringGeneration.dxy = 50;% dxy during generation, nm !!!!!!! fixed
PSFParameters.duringGeneration.xysize = round(PSFParameters.dxy*PSFParameters.xysize/PSFParameters.duringGeneration.dxy);%xysize during generation, pixels
PSFParameters.duringGeneration.centrePixelIndex = round(PSFParameters.duringGeneration.xysize/2);
if PSFParameters.duringGeneration.dxy>PSFParameters.lambda/4/PSFParameters.numAper
    warning('dxy SHOULD satisfy Nyquist sampling theorem when generating PSFs...');
end
PSFParameters.duringGeneration.dkxy = (2*pi)/(PSFParameters.duringGeneration.xysize*PSFParameters.duringGeneration.dxy);
PSFParameters.duringGeneration.kMax = (2*pi)/(2*(PSFParameters.lambda/2/PSFParameters.numAper));
PSFParameters.duringGeneration.kMaxPixels = PSFParameters.duringGeneration.kMax/PSFParameters.duringGeneration.dkxy;%radius
PSFParameters.defocus = (PSFParameters.defocusLower:PSFParameters.defocusInterval:PSFParameters.defocusUpper)*1000;%nm
% NA ratio
PSFParameters.lowNA = PSFParameters.numAper*PSFParameters.NARatio;
PSFParameters.duringGeneration.lowNAkMaxPixels = PSFParameters.duringGeneration.kMaxPixels * PSFParameters.NARatio;
% pinhole
PSFParameters.lowNADL = 920 / 2 / (1.05/10); %nm, ideal, fixed pinhole size in real world
% PSFParameters.lowNADL = PSFParameters.lambda / 2 / PSFParameters.lowNA; %nm, relatively fixed in different NARatio
PSFParameters.duringGeneration.apertureRadius = PSFParameters.lowNADL * PSFParameters.pinholeR_DL; %nm
PSFParameters.duringGeneration.apertureRadius_pixel = PSFParameters.duringGeneration.apertureRadius ./ PSFParameters.duringGeneration.dxy;
PSFParameters.duringGeneration.aperture = Circle([PSFParameters.duringGeneration.xysize,PSFParameters.duringGeneration.xysize],[PSFParameters.duringGeneration.centrePixelIndex,PSFParameters.duringGeneration.centrePixelIndex],[0,PSFParameters.duringGeneration.apertureRadius_pixel]);
PSFParameters.pinholeDefocus = 0;%nm
% angle distribution
PSFParameters.angDisty_pixel = round(PSFParameters.angDisty_normed' * PSFParameters.duringGeneration.kMaxPixels);
PSFParameters.angDistx_pixel = round(PSFParameters.angDistx_normed' * PSFParameters.duringGeneration.kMaxPixels);
PSFParameters.angDistGenerationMode = PSFParameters.angDistParameters.scanningMode;
PSFParameters.angDistx_normed = PSFParameters.angDistx_pixel / PSFParameters.duringGeneration.kMaxPixels;
PSFParameters.angDisty_normed = PSFParameters.angDisty_pixel / PSFParameters.duringGeneration.kMaxPixels;
% pupil phase
switch PSFParameters.pupilPhaseSet(1:2)
    case 'NO' % zero-phase
        PSFParameters.duringGeneration.pupilPlaneMatrix = zeros(PSFParameters.duringGeneration.xysize,PSFParameters.duringGeneration.xysize);
    case 'RA' % random phase
        stretch = str2double(PSFParameters.pupilPhaseSet(3:end));if isnan(stretch); stretch = 1;end
        polyMode = [4,6:15];
        polyCoef = randn(size(polyMode));
        PSFParameters.duringGeneration.pupilPlaneMatrix = zernikePhaseGenerator(PSFParameters.duringGeneration.xysize, PSFParameters.duringGeneration.kMaxPixels, savePath, polyMode, polyCoef, stretch);
    case 'SP' % spherical aberration only
        shpericalPhaseCoef = str2double(PSFParameters.pupilPhaseSet(3:end));if isnan(shpericalPhaseCoef); shpericalPhaseCoef = 1;end
        polyMode = [12];
        PSFParameters.duringGeneration.pupilPlaneMatrix = zernikePhaseGenerator(PSFParameters.duringGeneration.xysize, PSFParameters.duringGeneration.kMaxPixels, savePath, polyMode, shpericalPhaseCoef, 1);
    case 'SE' % estimated phase saved in folder ['systemParameters_path,'/pupilPhase/sePhase.mat']
        stretch = str2double(PSFParameters.pupilPhaseSet(3:end));if isnan(stretch); stretch = 1;end
        sePhase = loadFirstVariable([PSFParameters.systemParameters_path,'/pupilPhase/sePhase.mat']); % Estimated aber phase using DAO
        sePhase = imresize(stretch*sePhase,[ceil(2*PSFParameters.duringGeneration.kMaxPixels),ceil(2*PSFParameters.duringGeneration.kMaxPixels)],'bilinear');
        PSFParameters.duringGeneration.pupilPlaneMatrix = expandAndCentralizeAMatrix(sePhase,[PSFParameters.duringGeneration.xysize,PSFParameters.duringGeneration.xysize],0);
    otherwise %datestr(now, 'YYYYmmDD_HHMMSS'), use input phase
        if size(PSFParameters.duringGeneration.pupilPlaneMatrix,1) < PSFParameters.duringGeneration.xysize
            PSFParameters.duringGeneration.pupilPlaneMatrix = expandAndCentralizeAMatrix(PSFParameters.duringGeneration.pupilPlaneMatrix,[PSFParameters.duringGeneration.xysize,PSFParameters.duringGeneration.xysize],0);
        end
end
% K
PSFParameters.shift_K = PSFParameters.shiftKBias*1000 + [0, linspace(-PSFParameters.shiftKRange*1000/2,PSFParameters.shiftKRange*1000/2,PSFParameters.angleNum-1)];

%% %%%%%%%%%%%%%%%%%%%% generate new PSFs and save %%%%%%%%%%%%%%%%%%%%%%
% pupil to save at the end of this function
pupil_sum = Circle([PSFParameters.duringGeneration.xysize,PSFParameters.duringGeneration.xysize],[(PSFParameters.duringGeneration.xysize+1)/2,(PSFParameters.duringGeneration.xysize+1)/2],[0,PSFParameters.duringGeneration.kMaxPixels]);
pupil_ap_sum = Circle([PSFParameters.duringGeneration.xysize,PSFParameters.duringGeneration.xysize],[(PSFParameters.duringGeneration.xysize+1)/2,(PSFParameters.duringGeneration.xysize+1)/2],[0,PSFParameters.duringGeneration.kMaxPixels]);

% Get PSF one by one
disp('Initializing GPU...');
gpuDevice;
disp('GPU initialization finished');
disp(['caculate PSF start ||| Time : ',datestr(now, 'YYYYmmDD_HHMMSS')]);
for angleCount=1:PSFParameters.angleNum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [psfEx,pupil,pupil_lowNA,pupil_ap] = PSFSimuGPU(PSFParameters.duringGeneration.aperture,...
        PSFParameters.angDisty_pixel(angleCount),...
        PSFParameters.angDistx_pixel(angleCount),...
        PSFParameters.shift_K(angleCount),...
        PSFParameters.lambda,...
        PSFParameters.numAper,...
        PSFParameters.rindexObj,...
        PSFParameters.duringGeneration.dxy,...
        PSFParameters.defocus,...
        PSFParameters.duringGeneration.xysize,...
        PSFParameters.duringGeneration.dkxy,...
        PSFParameters.duringGeneration.kMaxPixels,...
        PSFParameters.duringGeneration.lowNAkMaxPixels,...
        PSFParameters.pinholeDefocus,...
        PSFParameters.duringGeneration.pupilPlaneMatrix,...
        PSFParameters.rindexSp,...
        PSFParameters.xysize);
    % two-photon excitation
    psf = abs(psfEx); psf(isnan(psf)) = 0;
    psf = abs(psf).^2; psf(isnan(psf)) = 0;
    % energy set
    switch PSFParameters.energySet
        case 0   % normalized with central slice 
            if angleCount == 1
                psfCentralSlice = psf(:,:,(end+1)/2);
                norEnergy = sum(psfCentralSlice(:));
            end
            psf = psf ./ norEnergy;
        case 1   % normalization for each layer
            psf=psf./sum(sum(psf,1),2);
            disp(['PSF normalized for each slice.']);
        case 2   % normalization for each angle
            psf=psf./sum(psf(:));
            disp(['PSF normalized for each angle.']);
        case 3   % normalization for each angle with their central slices
            psfCentralSlice = psf(:,:,(end+1)/2);
            psf=psf./sum(psfCentralSlice(:));
            disp(['PSF normalized for each angle with their central slices.']);
        case 4   % normalized with central slice, then set a minimum energy
            psfCentralSlice = psf(:,:,(end+1)/2);
            psf=psf./sum(psfCentralSlice(:));% normalization for this stack
            psfLayersEnergy = sum(sum(psf,1),2);
            max_psfLayersEnergy = max(psfLayersEnergy(:));
            minimumLayerEnergy = max_psfLayersEnergy/2;
            for zCount = 1:size(psf,3)
                psfThisLayer = psf(:,:,zCount);
                psfThisLayerEnergy = sum(psfThisLayer(:));
                if psfThisLayerEnergy < minimumLayerEnergy
                    psf(:,:,zCount) = psf(:,:,zCount)*minimumLayerEnergy/psfThisLayerEnergy;
                end
            end
            disp(['PSF normalized with their central slices, then set a minimum energy.']);
    end
    disp(['PSF ',num2str(angleCount),' generated. ||| Time : ',datestr(now, 'YYYYmmDD_HHMMSS')]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save every single psf as .tiff
    name_everyPSF=['psf_all_',num2str(angleCount)];
%     saveastiff_overwrite(psf,[savePath,'/',name_everyPSF,'.tiff']);
    psf_thisAngle = psf;
    save([savePath,'/',name_everyPSF,'.mat'],'psf_thisAngle');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% psf_sum & pupil of this psf
    psf_all(:,:,:,angleCount)=psf;
    pupil_sum = pupil_sum + imtranslate(pupil_lowNA,[PSFParameters.angDistx_pixel(angleCount),PSFParameters.angDisty_pixel(angleCount)]);
    pupil_ap_sum = pupil_ap_sum + abs(pupil_ap);
    PSFParameters.duringGeneration.pupil_ap_all(:,:,angleCount) = abs(pupil_ap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
PSFParameters.duringGeneration.pupil = pupil;

%% Save pupil and aperture(High&Low NA)
% saveastiff_overwrite(pupil_sum,[savePath '//pupil.tif']);
% saveastiff_overwrite(pupil_ap_sum,[savePath '//pupil_ap.tif']);
% saveastiff_overwrite(PSFParameters.duringGeneration.aperture,[savePath '//aperture.tif']);

%% Save psf_sum as .MAT for future convenience
disp('Saving psf_all as a .mat file...');
save([savePath,'//psf_all.mat'],'psf_all','-v7.3','-nocompression');
disp('Done!');

%% saving parameters during this PSF simulation process
save([savePath,'//PSFParameters.mat'],'PSFParameters','-v7.3');