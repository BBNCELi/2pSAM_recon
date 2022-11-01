% In-frame demotion algorithm for 2p-SAM
%ELi, 20220531
%ELi, 20220608, deconvRLGPU_fourier_DAO_forMotionCorr
function [proj_all, motionMap_sum] = demotion(psf_all,proj_all,PSFParameters,deconvOptions)
%% check inputs
[w,h,angleNum] = size(proj_all);
[w_,h_,s_,angleNum_] = size(psf_all);
if (w~=w_) || (h~=h_) || (angleNum~=angleNum_)
    error('Input error: psf_all and proj_all must have exactly the same size when performing in-frame motion correction...');
end
disp('In-frame motion correction start...');
%% preparations
motionCorrIter = 1;
thrld_meanEstimatedMotion_pixel = 0.01;
maxIter = 250;
updateWeight=ones(angleNum,1)*deconvOptions.equalWeight;
Xguess_init=ones(w_,h_,s_);
Xguess_init=Xguess_init./sum(Xguess_init(:))./angleNum;
motionMap = ones(angleNum, 2); MSE_iter = [];
motionMap_sum = zeros(angleNum, 2);
proj_all_ori = proj_all;
matSaveOn = 0; tifSaveOn = 0;
sidelobe = 100;
%% start
while (mean(motionMap(:).^2) > thrld_meanEstimatedMotion_pixel) && (motionCorrIter < maxIter) %pixel
    [Xguess_init_MCed,motionMap] = deconvRLGPU_fourier_DAO_forMotionCorr(psf_all,proj_all,PSFParameters,1,2,1,Xguess_init,updateWeight,deconvOptions.deconvSavepath,matSaveOn,tifSaveOn,sidelobe);
    motionMap_sum = motionMap_sum + motionMap;
    for angleNow = 1:angleNum
        proj_all(:,:,angleNow) = imtranslate(squeeze(proj_all_ori(:,:,angleNow)),[motionMap_sum(angleNow,2),motionMap_sum(angleNow,1)],'bilinear');
    end
    MSEHere = mean(motionMap(:).^2);
    MSE_iter(motionCorrIter) = MSEHere;
    disp(['Motion correction iter ',num2str(motionCorrIter),'... ', 'MSE = ',num2str(MSEHere),' pixels... thr_MSE: ',num2str(thrld_meanEstimatedMotion_pixel),' pixels...']);
%         saveastiff_overwrite(gather(Xguess_init_MCed),[deconvOptions.deconvSavepath,'//Xguess_prepMC_iter',num2str(motionCorrIter),'_MSE',num2str(mean(map_wavshape(:).^2)),'pixels.tiff']);
    motionCorrIter = motionCorrIter + 1;
end
%% save result
save([deconvOptions.deconvSavepath,'//inframe_motionCorrection_MSE_iter.mat'],'MSE_iter');
save([deconvOptions.deconvSavepath,'//inframe_motionCorrection_motionMap.mat'],'motionMap_sum');
saveastiff_overwrite(gather(Xguess_init_MCed),[deconvOptions.deconvSavepath,'//Xguess_MCed_iter',num2str(motionCorrIter),'_MSE',num2str(MSEHere),'pixels.tiff']);
end