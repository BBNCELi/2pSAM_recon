% 3D reconstruction with DAO using traditional RL algorithm
% 
% NOTES
%     The traditional RL algorithm for light-field reconstruction (also
%     termed as 3D deconvolution). We here implement it on GPU but compare
%     to fRL, the memory usage (both GPU and CPU memory) are relatively
%     high and the speed is 2-10 times slow.
% 
% REFERENCE
%     Broxton M, Grosenick L, Yang S, et al. Wave optics theory and 3-D
%     deconvolution for the light field microscope [J]. Optics Express, 
%     2013, 21(21): 25418-25439.
% 
% ELi, 20230131, rearrange inputs
function [Xguess,shiftMap,dispMap] = recon_RL_GPU(psfs,projs,reconOpts,Xguess)
% INPUT
%     psfs      - 3D psfs used for reconstrution, (usually) 4D matrix
%     projs     - 2D projections to be reconstructed, (usually) 3D matrix
%     reconOpts - options for 3D reconstruction, struct (default)
%                 .maxIter: max iterations (10)
%                 .upWeight: weight for volume update in each iteration (.5)
%                 .upSeq: updating sequence for projections from different angles
%                 .CAIndex: index of centre projection in projs (1)
%                 .savepath: folder to save results. If not given, generate a folder in saveFolder and save results there
%                 .tifSaveOn: iterations save flag (.tif), 0/1 (0)
%                 .matSaveOn: iterations save flag (.mat), 0/1 (0)
%                 .dispOn: display reconstruction process flag, 0/1 (1)
%                 .DAO: DAO flag, 0/1
%                 .DAOOpts: options for DAO
%                     .maxIter_init: max iteration without DAO
%                     .shiftEstFunc: function used for shift estimation
%                     .patchN: number of patches for multi-site shift estimation
%                     .patchOvFactor: overlap factor for multi-site shift estimation
%                     .minPatchSize: minimum patch size for multi-site shift estimation
%                     .sidelobe: sidelobe in pixel for shift estimation
%                     .sidelobez: use only the high-resolution regime for shift estimation
%     Xguess    - input volume

%% check inputs and preparations
%%%size
[psf_r,psf_c,psf_s,angleNum_] = size(psfs);
[proj_r,proj_c,angleNum] = size(projs);
if angleNum ~= angleNum_; error('INPUT ERROR: angles in psfs and projs should be the same'); end
if (psf_r > proj_r) || (psf_c > proj_c); error('INPUT ERROR: projection should be larger than psf'); end
%%%GPU
Xguess=gpuArray(single(Xguess)); %use single to speed up
dispOn = reconOpts.dispOn;
%%%psf_t i.e. the backward projection kernel in RL deconvolution
if dispOn; disp('Calculating psf_t...'); end
psf_t = zeros(psf_r,psf_c,psf_s,angleNum,'single');
for sliceCount=1:psf_s
    for psfCount=1:angleNum
        psf_t(:,:,sliceCount,psfCount)=rot90(squeeze(psfs(:,:,sliceCount,psfCount)),2);
    end
end
%%%htf
if dispOn; disp('Calculating htf...'); end
htf = gpuArray.zeros(proj_r,proj_c,psf_s,angleNum,'single');
for i = 1:angleNum
    htf(:,:,:,i) = backwardProj_RL_GPU(gpuArray(squeeze(psf_t(:,:,:,i))),gpuArray.ones(proj_r,proj_c,'single'));
end
%%%extract struct
maxIter = reconOpts.maxIter;
maxIter_init = maxIter;
upWeight = reconOpts.upWeight;
upSeq = reconOpts.upSeq;
CAIndex = reconOpts.CAIndex;
savepath = reconOpts.savepath;
tifSaveOn = reconOpts.tifSaveOn;
matSaveOn = reconOpts.matSaveOn;
dispOn = reconOpts.dispOn;
shiftMap = 0;
dispMap = 0;
%%%DAO
DAO = reconOpts.DAO;
if DAO
    if ~isfield(reconOpts,'DAOOpts');reconOpts.DAOOpts=struct;end
    DAOOpts = reconOpts.DAOOpts;
    if ~isfield(DAOOpts,'maxIter_init');DAOOpts.maxIter_init = 1;end
    maxIter_init = DAOOpts.maxIter_init;
    if ~isfield(DAOOpts,'shiftEstFunc');DAOOpts.shiftEstFunc='corr';end
    if ~isfield(DAOOpts,'patchN');DAOOpts.patchN=1;end
    if ~isfield(DAOOpts,'patchOvFactor');DAOOpts.patchOvFactor=0.1;end
    if ~isfield(DAOOpts,'minPatchSize');DAOOpts.minPatchSize=32;end
    if ~isfield(DAOOpts,'sidelobe');DAOOpts.sidelobe=20;end
    if ~isfield(DAOOpts,'sidelobez');DAOOpts.sidelobez=30;end
    shiftMap = zeros(DAOOpts.patchN,DAOOpts.patchN,2,angleNum);
    dispMap = zeros(proj_r, proj_c, 2, angleNum, 'single');
end

%% update volume with or without DAO
for iterNow = 1:maxIter
    %%% DAO
    if DAO == 1 && iterNow > maxIter_init
        if dispOn
            fprintf('Estimating shift and disparity map...');
        end
        for angleNow = 1:angleNum
            %%%prepare projection in this angle
            projNow = gpuArray(squeeze(projs(:,:,angleNow)));

            %%%forward project to get HXguess
            XguessDAO = Xguess;
            XguessDAO(:,:,[1:DAOOpts.sidelobez,psf_s-DAOOpts.sidelobez:psf_s])=0;
            HXguess = forwardProj_RL_GPU(gpuArray(squeeze(psfs(:,:,:,angleNow))),XguessDAO);
            clear XguessDAO

            %%%estimate shift and disparity map
            [shiftMap(:,:,:,angleNow),dispMap(:,:,:,angleNow)] = ...
                shiftEst(projNow, HXguess, DAOOpts);
        end
        %%set centre angle shift to 0
        dispMap = dispMap - dispMap(:,:,:,CAIndex);
        shiftMap = shiftMap - shiftMap(:,:,:,CAIndex);
        if dispOn
            fprintf(' done\n');
        end
    end

    %%% reconstruction
    for upCount = 1:length(upSeq)
        upAngles = upSeq{upCount};
        upAnglesNum = length(upAngles);

        for angleNow = upAngles
            %%%prepare projection in this angle
            projNow = squeeze(projs(:,:,angleNow));
            if DAO == 1 && iterNow > maxIter_init
                projNow = imtranslate_disparityMap(projNow,dispMap(:,:,:,angleNow));
            end
            projNow = gpuArray(projNow);

            %%%forward project to get HXguess
            HXguess = forwardProj_RL_GPU(gpuArray(squeeze(psfs(:,:,:,angleNow))),Xguess);

            %%%calculate error and backward project
            errorEM = projNow./HXguess;errorEM(~isfinite(errorEM))=0;
            if ~exist('errorBack','var'); errorBack = gpuArray.zeros(proj_r,proj_c,psf_s,'single'); end
            errorBack = backwardProj_RL_GPU(gpuArray(squeeze(psf_t(:,:,:,angleNow))),errorEM) ./ htf(:,:,:,angleNow) * 1/upAnglesNum + errorBack;
        end

        %%%update Volume
        Xguess = Xguess.*errorBack.*upWeight+(1-upWeight).*Xguess;
            
        %%%finish
        clear errorBack;%To avoid out of GPU memory
        Xguess(isnan(Xguess)) = 0; Xguess = real(Xguess); Xguess(Xguess<0) = 0;
        if dispOn
            disp(['Reconstructing volume... ',...
                array2str(upAngles,'+'),'th ||| ', num2str(angleNum),' angle...',...
                num2str(iterNow),'th ||| ', num2str(maxIter),' iter...',...
                ' Energy=',num2str(sum(Xguess(:))), '...',...
                ' TimeNow: ',datestr(now,'YYYYmmDD_HHMMSS')]);
        end
    end

    %%% Save
    stackTempForSave = gather(Xguess);
    saveastiff_overwrite(stackTempForSave, [savepath,'/Xguess_iter',num2str(iterNow),'.tif'],matSaveOn,tifSaveOn);
end