% Forward-project a volume with known PSFs

%%ELiiiiiii, 20211211
function proj_all = forwardProjectAndSave(sample,psf_all,savepath,speed_pixels)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w,h,s] = size(sample);
[w_,h_,s_,angleNum] = size(psf_all);
if (w~=w_) || (h~=h_) || (s~=s_)
    error('Input error: sample and PSFs must have exactly the same size');
end
if ~exist('speed_pixels','var'); speed_pixels = 0; end
%% %%%%%%%%%%%%%%%%%%%%%%% convolve slice by slice %%%%%%%%%%%%%%%%%%%%%%%
proj_all = zeros(w, h, angleNum, 'single');
for psfCount = 1:angleNum
    disp(['Imaging from the ', num2str(psfCount),'th angle || Timenow: ',datestr(now, 'YYYYmmDD_HHMMSS')]);
    psfsingle=psf_all(:,:,:,psfCount);
    % moving sample
    if (psfCount ~= 1) && (speed_pixels ~= 0)
        sample = imtranslate(sample,[0,speed_pixels,0],'bilinear');
    end

    %%%%%%%%%%
    %%% method1: conv2 ~180s
%     projection = zeros(w, h, 'single');
%     for sCount = 1:s
%         projection = projection + conv2(sample(:,:,sCount), psfsingle(:,:,sCount), 'same');
%     end
    %%% method2: fft ~1s
    projection = zeros(w, h, 'single');
    for sCount = 1:s
        projection = projection + real(ifftForImage(fftForImage(sample(:,:,sCount)) .* fftForImage(psfsingle(:,:,sCount))));
    end
    %%% method3: fft-GPU ~0.5s
%     projection = gpuArray.zeros(w, h, 'single');
%     psfsingle = gpuArray(psfsingle);
%     sample = gpuArray(sample);
%     for sCount = 1:s
%         projection = projection + real(ifftForImage(fftForImage(sample(:,:,sCount)) .* fftForImage(psfsingle(:,:,sCount))));
%     end
%     projection = gather(projection);
    %%%%%%%%%%

    proj_all(:,:,psfCount) = projection;
end
%% %%%%%%%%%%%%%%%%%%%%%%% save if savepath is given %%%%%%%%%%%%%%%%%%%%%%%
if exist('savepath','var')
    saveastiff_overwrite(proj_all, [savepath,'//proj_all.tif']);
    save([savepath,'//proj_all.mat'],'proj_all','-mat','-v7.3');
end