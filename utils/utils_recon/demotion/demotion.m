% In-frame motion correction for 2p-SAM
%
% ELi, 20230131, rearrange inputs
function [projs, motionMap, dispMap] = demotion(psfs,projs,reconOpts)
% INPUT
%     psfs      - 3D psfs used for motion correction, (usually) 4D matrix
%     projs     - 2D projections to be motion correction, (usually) 3D matrix
%     reconOpts - options for 3D reconstruction, struct (default)
%                 .solver: reconstruction algorithm (fRL)
%                 .upWeight: weight for volume update in each iteration (.5)
%                 .upSeq: updating sequence for projections from different angles
%                 .CAIndex: index of centre projection in projs (1)
%                 .savepath: folder to save results. If not given, generate a folder in saveFolder and save results there
%                 .demotion: demotion flag, 0/1
%                 .demotionOpts: options for motion correction
%                     .maxIter_init: iterations without shift estimation in each iteration of demotion
%                     .maxIter: iterations without+out shift estimation in each iteration of demotion
%                     .shiftEstFunc: function used for shift estimation
%                     .patchN: number of patches for multi-site shift estimation
%                     .patchOvFactor: overlap factor for multi-site shift estimation
%                     .minPatchSize: minimum patch size for multi-site shift estimation
%                     .sidelobe: sidelobe in pixel for shift estimation
%                     .sidelobez: use only the high-resolution regime for shift estimation
%                     .thrMSE: the motion map converges if its increments are small enough
%                     .maxIterDemotion: max iterations to quit demotion if the motion map does not converge
%                     .maxShift: max shift in pixel for each iteration of demotion
%                     .angDistx_normed: PSF angle distribution, x, to remove the defocus item when initializing motion map
%                     .angDisty_normed: PSF angle distribution, y, to remove the defocus item when initializing motion map

%% check inputs and preparations
%%%size
[proj_r,proj_c,angleNum] = size(projs);
[psf_r,psf_c,psf_s,angleNum_] = size(psfs);
if (proj_r~=psf_r) || (proj_c~=psf_c) || (angleNum~=angleNum_)
    error('Input error: psfs and projs must have exactly the same size when performing in-frame motion correction...');
end
%%%default options
if ~isfield(reconOpts,'demotionOpts');reconOpts.demotionOpts=struct;end
if ~isfield(reconOpts.demotionOpts,'maxIter_init');reconOpts.demotionOpts.maxIter_init = 3;end
if ~isfield(reconOpts.demotionOpts,'maxIter');reconOpts.demotionOpts.maxIter = 4;end
if ~isfield(reconOpts.demotionOpts,'shiftEstFunc');reconOpts.demotionOpts.shiftEstFunc = 'imregtform';end
if ~isfield(reconOpts.demotionOpts,'patchN');reconOpts.demotionOpts.patchN = 1;end
if ~isfield(reconOpts.demotionOpts,'patchOvFactor');reconOpts.demotionOpts.patchOvFactor = 0.1;end
if ~isfield(reconOpts.demotionOpts,'minPatchSize');reconOpts.demotionOpts.minPatchSize = 32;end
if ~isfield(reconOpts.demotionOpts,'sidelobe');reconOpts.demotionOpts.sidelobe = 100;end
if ~isfield(reconOpts.demotionOpts,'sidelobez');reconOpts.demotionOpts.sidelobez = 35;end
if ~isfield(reconOpts.demotionOpts,'thrMSE');reconOpts.demotionOpts.thrMSE = 0.01;end
if ~isfield(reconOpts.demotionOpts,'maxIterDemotion');reconOpts.demotionOpts.maxIterDemotion = 100;end
if ~isfield(reconOpts.demotionOpts,'maxShift');reconOpts.demotionOpts.maxShift = inf;end
%%%Copy demotion options to DAO options so that we don't have to write a
%%%demotion reconstruction function particularly
reconOpts.maxIter = reconOpts.demotionOpts.maxIter;
reconOpts.DAO = 1;
reconOpts.DAOOpts = reconOpts.demotionOpts;
reconOpts.tifSaveOn = 0;
reconOpts.matSaveOn = 0;
reconOpts.dispOn = 0;
%%%extract struct
thrMSE = reconOpts.demotionOpts.thrMSE;
maxIter = reconOpts.demotionOpts.maxIterDemotion;
maxShift = reconOpts.demotionOpts.maxShift;
angDistx_normed = reconOpts.demotionOpts.angDistx_normed;
angDisty_normed = reconOpts.demotionOpts.angDisty_normed;
%%%to be saved
MSE_iter = [];
motionMap = zeros(reconOpts.demotionOpts.patchN,reconOpts.demotionOpts.patchN,2,angleNum);%motion map
dispMap = zeros(proj_r, proj_c, 2, angleNum, 'single');%disparity map

%%%save a copy of projs
projs_ori = projs; %original projections
%%%volume initialization
Xguess_init = initXguess(psfs,projs);

disp('In-frame motion correction start...');

%% initialization for large motions
% Before the iterations of motion correction, we add an initialization
% process for motionMap and dispMap.
% 
% We estimate the shifts between projections and the projection at the
% centre view (CA) to get a shift map. This shift map contains 3 
% components: shifts induced by motion (which is what we want), shifts
% induced by defocus (which is the core information for volume
% reconstruction) and shifts induced by aberrations (which could be
% estimated later by multi-site DAO methods; note that while tilting and
% defocus can also be viewed as low-order aberrations, when talking about
% aberrations here we mean the higher-order ones). When shifts induced by
% motion and defocus are comparativly large (tens of pixels between
% projections), aberration-induced shifts are usually small (several pixels
% between projections) and inhomogenous and thus can be ignored. At the
% same time, as the shifts induced by defocus have a fixed distribution at
% different viewpoints and the viewpoints have already been given when
% modeling PSFs, we can directly (and roughly) estimate the defocus-induced
% shifts and then subtract it from the shift map. We use this
% defocus-removed shift map to initialize the motion map and the disparity
% map. This initializaiton greatly speeds up the convergence especially
% when large motions are present.

largeMotions = 1;
if largeMotions
    motionMap_init = shiftToCA(projs_ori,reconOpts.CAIndex,'corr'); %we use the simplest correlation 
    motionMap_init = removeDefocusItem(motionMap_init,angDistx_normed,angDisty_normed,reconOpts.CAIndex);
    motionMap(:,:,1,:) = motionMap_init(1,:);
    motionMap(:,:,2,:) = motionMap_init(2,:);
    for angleNow = 1:angleNum
        dispMap(:,:,1,angleNow) = motionMap(:,:,1,angleNow);
        dispMap(:,:,2,angleNow) = motionMap(:,:,2,angleNow);
        projs(:,:,angleNow) = imtranslate(squeeze(projs_ori(:,:,angleNow)),[motionMap_init(2,angleNow),motionMap_init(1,angleNow)],'bilinear');
    end
end

%% demotion iterations
iter = 1;
while (iter <= maxIter) %stop at maxIter if not converge
    %%%estimate motion in reconstruction
    [~,motionMap_incre,dispMap_incre] = recon_fRL_GPU(psfs,projs,reconOpts,Xguess_init);

    %%%max motion in each iteration
    motionMap_incre(motionMap_incre > maxShift) = maxShift;
    motionMap_incre(motionMap_incre < -maxShift) = -maxShift;
    motionMap = motionMap + motionMap_incre;
    dispMap_incre(dispMap_incre > maxShift) = maxShift;
    dispMap_incre(dispMap_incre < -maxShift) = -maxShift;
    dispMap = dispMap + dispMap_incre;

    %%%move projs according to disparity map now
    for angleNow = 1:angleNum
        projs(:,:,angleNow) = imtranslate_disparityMap(squeeze(projs_ori(:,:,angleNow)),dispMap(:,:,:,angleNow));
    end

    %%%converge or not
    MSEHere = mean(motionMap_incre(:).^2);
    MSE_iter(iter) = MSEHere;
    disp(['Motion correction iter ',num2str(iter),'... ', 'MSE = ',num2str(MSEHere),' ... thr_MSE: ',num2str(thrMSE),' ...']);
    iter = iter + 1;
    if MSEHere < thrMSE
        break;
    end
end
disp('In-frame motion correction finished');

%% save result
save([reconOpts.savepath,'//demotion',datestr(now, 'YYYYmmDD_HHMMSS'),'.mat'],'motionMap','dispMap','MSE_iter','reconOpts');
end