% Fourier RL deconvolution with DAO, using GPU
%ELi, 20220531, final
function [Xguess,shifts_estimated] = deconvRLGPU_fourier_DAO(psf_all,proj_all,PSFParameters,maxIter,maxIter_init,AO,Xguess,updateWeight,savepath,matSaveOn,tifSaveOn,sidelobe)
%% check inputs and preparations
[psf_r,psf_c,psf_s,angleNum_] = size(psf_all);
[img_r,img_c,angleNum] = size(proj_all);
if angleNum ~= angleNum_
    error('Mismatch between PSF and blur_image size!');
end
disp('Initializing GPU...');
Xguess=gpuArray(single(Xguess)); largepsf = gpuArray.zeros(img_r,img_c,psf_s,'single');
disp('GPU initialization finished');
%% iteration sequence settings
% Angle-updating sequence in iterations
[psf_iterSequence,centreAngleIndex] = getIterSequence(PSFParameters);
% proj_all(:,:,centreAngleIndex) is the centre angle image
% psf_all(:,:,:,centreAngleIndex) is the centre angle psf

% Angle-updating sequence when initializing volume
initializeXguessAnglesIndex = psf_iterSequence;%someOtherOtherfunction(psfNumberInTotal);

%% Initialization
for iterNow = 1:maxIter_init
    for i = 1:length(initializeXguessAnglesIndex)
        %%%preparations
        largepsf((img_r-psf_r)/2+1:(img_r+psf_r)/2,(img_c-psf_c)/2+1:(img_c+psf_c)/2,:) = gpuArray(squeeze(psf_all(:,:,:,initializeXguessAnglesIndex(i))));
        projHere = squeeze(proj_all(:,:,initializeXguessAnglesIndex(i)));
        
        %%%Forward and backward projection
        HXguessFFT = forwardProjectRLGPU_fourier(largepsf,Xguess,psf_s);
        errorBack = backwardProjectRLGPU_fourier(largepsf,HXguessFFT,psf_s,img_r,img_c,projHere);
    
        %%%Update Volume
        Xguess = Xguess.*errorBack.*updateWeight(initializeXguessAnglesIndex(i))+(1-updateWeight(initializeXguessAnglesIndex(i))).*Xguess;
        
        %%%finish
        clear errorBack;%To save GPU memory
        Xguess(isnan(Xguess)) = 0;Xguess = real(Xguess);Xguess(Xguess<0) = 0;
%         disp(['Initializing Xguess... ',...
%             num2str(initializeXguessAnglesIndex(i)),'th ||| ', num2str(angleNum),' angle...',...
%             num2str(iterNow),'th ||| ', num2str(maxIter_init),' iter...',...
%             ' Energy=',num2str(sum(Xguess(:))), '...',...
%             ' TimeNow: ',datestr(now,'YYYYmmDD_HHMMSS')]);
    end
end
    disp(['Xguess initialized... ',...
        ' Energy=',num2str(sum(Xguess(:))), '...',...
        ' TimeNow: ',datestr(now,'YYYYmmDD_HHMMSS')]);
    %%%Save
    stackTempForSave = gather(Xguess);
    saveastiff_overwrite(stackTempForSave, [savepath,'//Xguess_iter0.tif'],matSaveOn,tifSaveOn);

%% iters
%%%DAO Preparations
shifts_estimated = zeros(angleNum, 2);
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
            %%%preparations
            largepsf((img_r-psf_r)/2+1:(img_r+psf_r)/2,(img_c-psf_c)/2+1:(img_c+psf_c)/2,:)=gpuArray(squeeze(psf_all(:,:,:,angleNow)));
            projHere = imtranslate(squeeze(proj_all(:,:,angleNow)),[shifts_estimated(angleNow,2),shifts_estimated(angleNow,1)],'bilinear');
            
            %%%Forward and backward projection
            HXguessFFT = forwardProjectRLGPU_fourier(largepsf,Xguess,psf_s);
            errorBack = backwardProjectRLGPU_fourier(largepsf,HXguessFFT,psf_s,img_r,img_c,projHere);

            %%Update volume
            Xguess = Xguess.*errorBack.*updateWeight(angleNow)+(1-updateWeight(angleNow)).*Xguess;
            
            %%%Finish
            clear errorBack;%To save GPU memory
            Xguess(isnan(Xguess)) = 0;Xguess = real(Xguess);Xguess(Xguess<0) = 0;
%             disp(['Updating Xguess... ',...
%                 num2str(angleNow),'th ||| ', num2str(angleNum),' angle...',...
%                 num2str(iterNow),'th ||| ', num2str(maxIter),' iter...',...
%                 ' Energy=',num2str(sum(Xguess(:))), '...',...
%                 ' TimeNow: ',datestr(now,'YYYYmmDD_HHMMSS')]);
        end
    end
    disp(['Xguess updated... ',...
        num2str(iterNow),'th ||| ', num2str(maxIter),' iter...',...
        ' Energy=',num2str(sum(Xguess(:))), '...',...
        ' TimeNow: ',datestr(now,'YYYYmmDD_HHMMSS')]);

    %% DAO
    if AO == 1 || AO == 2 %Do DAO only if AOStar == 1 or 2
%         disp('Estimating shifts...');
        for angleNow = psf_iterSequence
            %%%preparations
            largepsf((img_r-psf_r)/2+1:(img_r+psf_r)/2,(img_c-psf_c)/2+1:(img_c+psf_c)/2,:)=gpuArray(squeeze(psf_all(:,:,:,angleNow)));
            
            %%%Forward projection and get HXguess
            HXguessFFT = forwardProjectRLGPU_fourier(largepsf,Xguess,psf_s);
            HXguess = abs(ifftForImage(HXguessFFT));
            
            %%%DAO
            sub_HXguess=HXguess(yRangeForCorr,xRangeForCorr);
            sub_proj=gpuArray(squeeze(proj_all(yRangeForCorr,xRangeForCorr,angleNow)));
%             [shifts_estimated(angleNow,1),shifts_estimated(angleNow,2)] = findImgShiftByPhaseCorr(sub_proj,sub_HXguess);
            [shifts_estimated(angleNow,1),shifts_estimated(angleNow,2)] = findImgShift(sub_proj,sub_HXguess);
%             [shifts_estimated(angleNow,1),shifts_estimated(angleNow,2)] = findImgShift_imregtform(sub_proj,sub_HXguess);
        end

        %%% shifts_estimated modification
        shifts_estimated_centreAngleIndex_xDirection=shifts_estimated(centreAngleIndex,1);
        shifts_estimated_centreAngleIndex_yDirection=shifts_estimated(centreAngleIndex,2);
        shifts_estimated(:,1)=(squeeze(shifts_estimated(:,1))-shifts_estimated_centreAngleIndex_xDirection).*mask;
        shifts_estimated(:,2)=(squeeze(shifts_estimated(:,2))-shifts_estimated_centreAngleIndex_yDirection).*mask;

        %%%defocus
        if AO == 2
            k1 = PSFParameters.angDistx_normed.*squeeze(shifts_estimated(:,2)).*mask+PSFParameters.angDisty_normed.*squeeze(shifts_estimated(:,1)).*mask;
            k2 = PSFParameters.angDisty_normed.*PSFParameters.angDisty_normed.*mask+PSFParameters.angDistx_normed.*PSFParameters.angDistx_normed.*mask;
            k=sum(k1(:))/sum(k2(:));
            shifts_estimated(:,1)=squeeze(shifts_estimated(:,1))+k*PSFParameters.angDisty_normed;
            shifts_estimated(:,2)=squeeze(shifts_estimated(:,2))+k*PSFParameters.angDistx_normed;
            shifts_estimated_centreAngleIndex_xDirection=shifts_estimated(centreAngleIndex,1);
            shifts_estimated_centreAngleIndex_yDirection=shifts_estimated(centreAngleIndex,2);
            shifts_estimated(:,1)=(squeeze(shifts_estimated(:,1))-shifts_estimated_centreAngleIndex_xDirection).*mask;
            shifts_estimated(:,2)=(squeeze(shifts_estimated(:,2))-shifts_estimated_centreAngleIndex_yDirection).*mask;
        end
%         save([savepath,'//shift_estimated_iter_',num2str(iterNow),'.mat'],'shifts_estimated');
    end

    %% Save
    stackTempForSave = gather(Xguess);
    saveastiff_overwrite(stackTempForSave, [savepath,'/Xguess_iter',num2str(iterNow),'.tif'],matSaveOn,tifSaveOn);
end