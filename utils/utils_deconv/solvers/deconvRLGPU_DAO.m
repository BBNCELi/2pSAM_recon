% RL deconvolution with DAO, using GPU
%ELi, 20220531, final
function [Xguess,shifts_estimated] = deconvRLGPU_DAO(psf_all,proj_all,PSFParameters,maxIter,maxIter_init,AO,Xguess,updateWeight,savepath,matSaveOn,tifSaveOn,sidelobe)
[psf_r,psf_c,psf_s,angleNum_] = size(psf_all);
[img_r,img_c,angleNum] = size(proj_all);
if angleNum ~= angleNum_
    error('Mismatch between PSF and blur_image size!');
end
Xguess=gpuArray(single(Xguess));
%Rotate 180 degrees for each slice of psf to get psf_t, i.e. the kernel for RL backward projection.
disp('Calculating psf_t...');
psf_t = zeros(psf_r,psf_c,psf_s,angleNum,'single');
for sliceCount=1:psf_s
    for psfCount=1:angleNum
        psf_t(:,:,sliceCount,psfCount)=rot90(squeeze(psf_all(:,:,sliceCount,psfCount)),2);
    end
end
%get htf
disp('Calculating htf...');
htf = gpuArray.zeros(img_r,img_c,psf_s,angleNum,'single');
for i = 1:angleNum
    htf(:,:,:,i) = backwardProjectGPU(gpuArray(squeeze(psf_t(:,:,:,i))),gpuArray.ones(img_r,img_c,'single'));
end

%% iteration sequence settings
% Angle-updating sequence in iterations
[psf_iterSequence,centreAngleIndex] = getIterSequence(PSFParameters);
% proj_all(:,:,centreAngleIndex) is the centre angle image
% psf_all(:,:,:,centreAngleIndex) is the centre angle psf

% Angle-updating sequence when initializing volume
initializeXguessAnglesIndex = psf_iterSequence;%someOtherOtherfunction(psfNumberInTotal);


%% Initializations
for iterNow = 1:maxIter_init
    for i = 1:length(initializeXguessAnglesIndex)
        %%%preparations
        projHere = squeeze(proj_all(:,:,initializeXguessAnglesIndex(i)));
        
        %%%Forward and backward projection
        HXguess = forwardProjectGPU(gpuArray(squeeze(psf_all(:,:,:,initializeXguessAnglesIndex(i)))),Xguess);
        errorEM = projHere./HXguess; errorEM(~isfinite(errorEM))=0;
        XguessCor = backwardProjectGPU(gpuArray(squeeze(psf_t(:,:,:,initializeXguessAnglesIndex(i)))),errorEM) ./ htf(:,:,:,initializeXguessAnglesIndex(i));
        
        %%%Update volume
        Xguess=Xguess.*XguessCor.*updateWeight(initializeXguessAnglesIndex(i))+(1-updateWeight(initializeXguessAnglesIndex(i))).*Xguess;
        
        %%%Finish
        Xguess(isnan(Xguess)) = 0;Xguess(Xguess<0) = 0;
        disp(['Initializing Xguess... ',...
            num2str(initializeXguessAnglesIndex(i)),'th ||| ', num2str(angleNum),' angle...',...
            num2str(iterNow),'th ||| ', num2str(maxIter_init),' iter...',...
            ' Energy=',num2str(sum(Xguess(:))), '...',...
            ' TimeNow: ',datestr(now,'YYYYmmDD_HHMMSS')]);
    end
end    
    %%%Save
    stackTempForSave = gather(Xguess);
    saveastiff_overwrite(stackTempForSave, [savepath,'//XguessRL_iter0.tif'], matSaveOn, tifSaveOn)

%% iters
%%%DAO Preparations
shifts_estimated=zeros(angleNum,2);
if ~exist('sidelobe','var'); sidelobe = 20; end
yRangeForCorr = sidelobe + 1 : img_r - sidelobe; xRangeForCorr = sidelobe + 1 : img_c - sidelobe;
mask=ones(angleNum,1);

%%%iters start
for iterNow = 1:maxIter
    %% Deconv
    for angleNow = psf_iterSequence
        if updateWeight(angleNow)==0
            disp(['  iter ' num2str(iterNow) ' ||| ' num2str(maxIter),' (angleNow=',num2str(angleNow),'  continue ....' ]);
            continue;
        else
            %%%Preparations
            projHere = imtranslate(squeeze(proj_all(:,:,angleNow)),[shifts_estimated(angleNow,2),shifts_estimated(angleNow,1)],'cubic');
            
            %%%Forward and backward projection
            HXguess = forwardProjectGPU(gpuArray(squeeze(psf_all(:,:,:,angleNow))),Xguess);
            errorEM = projHere./HXguess;errorEM(~isfinite(errorEM))=0;
            XguessCor = backwardProjectGPU(gpuArray(squeeze(psf_t(:,:,:,angleNow))),errorEM) ./ htf(:,:,:,angleNow);
            
            %%%Update volume
            Xguess=Xguess.*XguessCor.*updateWeight(angleNow)+(1-updateWeight(angleNow)).*Xguess;
            
            %%%Finish
            Xguess(isnan(Xguess)) = 0;Xguess(Xguess<0) = 0;
            disp(['Updating Xguess... ',...
                num2str(angleNow),'th ||| ', num2str(angleNum),' angle...',...
                num2str(iterNow),'th ||| ', num2str(maxIter),' iter...',...
                ' Energy=',num2str(sum(Xguess(:))), '...',...
                ' TimeNow: ',datestr(now,'YYYYmmDD_HHMMSS')]);
        end
    end

    %% DAO
    if AO == 1 || AO == 2 %Do DAO only if AOStar == 1 or 2
        for angleNow = psf_iterSequence
            %%%Forward projection
            HXguess = forwardProjectGPU(gpuArray(squeeze(psf_all(:,:,:,angleNow))),Xguess);
            
            %%%DAO
            sub_HXguess=HXguess(yRangeForCorr,xRangeForCorr);
            sub_proj=gpuArray(squeeze(proj_all(yRangeForCorr,xRangeForCorr,angleNow))); 
%             [shifts_estimated(angleNow,1),shifts_estimated(angleNow,2)] = findImgShiftByPhaseCorr(sub_proj,sub_HXguess);
            [shifts_estimated(angleNow,1),shifts_estimated(angleNow,2)] = findImgShift(sub_proj,sub_HXguess);
%             [shifts_estimated(angleNow,1),shifts_estimated(angleNow,2)] = findImgShift_imregtform(sub_proj,sub_HXguess);
        end
        %%% map_wavshape modification
        shifts_estimated_centreAngleIndex_xDirection=shifts_estimated(centreAngleIndex,1);
        shifts_estimated_centreAngleIndex_yDirection=shifts_estimated(centreAngleIndex,2);
        shifts_estimated(:,1)=(squeeze(shifts_estimated(:,1))-shifts_estimated_centreAngleIndex_xDirection).*mask;
        shifts_estimated(:,2)=(squeeze(shifts_estimated(:,2))-shifts_estimated_centreAngleIndex_yDirection).*mask;

        %%%defocus
        if AO == 2
            k1 = PSFParameters.shift_x_normalized.*squeeze(shifts_estimated(:,2)).*mask+PSFParameters.shift_y_normalized.*squeeze(shifts_estimated(:,1)).*mask;
            k2 = PSFParameters.shift_y_normalized.*PSFParameters.shift_y_normalized.*mask+PSFParameters.shift_x_normalized.*PSFParameters.shift_x_normalized.*mask;
            k=sum(k1(:))/sum(k2(:));
            shifts_estimated(:,1)=squeeze(shifts_estimated(:,1))+k*PSFParameters.shift_y_normalized;
            shifts_estimated(:,2)=squeeze(shifts_estimated(:,2))+k*PSFParameters.shift_x_normalized;
            shifts_estimated_centreAngleIndex_xDirection=shifts_estimated(centreAngleIndex,1);
            shifts_estimated_centreAngleIndex_yDirection=shifts_estimated(centreAngleIndex,2);
            shifts_estimated(:,1)=(squeeze(shifts_estimated(:,1))-shifts_estimated_centreAngleIndex_xDirection).*mask;
            shifts_estimated(:,2)=(squeeze(shifts_estimated(:,2))-shifts_estimated_centreAngleIndex_yDirection).*mask;
        end
%         save([savepath,'/shift_estimated_iter_',num2str(iterNow),'.mat'],'shifts_estimated');
    end
    
    %% Save
    stackTempForSave = gather(Xguess);
    saveastiff_overwrite(stackTempForSave, [savepath,'//XguessRL_iter',num2str(iterNow),'.tif'], matSaveOn, tifSaveOn);
end