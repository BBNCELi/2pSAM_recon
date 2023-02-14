% Backward projector in fRL for 3D reconstruction
% 
% REFERENCE
%     Wu J, et al. Iterative tomography with digital adaptive optics permits
%     hour-long intravital observation of 3D subcellular dynamics at millisecond
%     scale [J]. Cell, 2021, 184(12): 3318-3332 e3317.
% 
% ELi, 20230131, add comments
function errorBack = backwardProj_fRL_GPU(largepsf,HXguessFFT,projNow)
% INPUT
%     largepsf   - psf but in the size of projection
%     HXguessFFT - Fourier transform of HXguess, given by forwardProj_fRL_GPU.m
%     projNow    - projection to be used

%% size
[proj_r,proj_c,psf_s] = size(largepsf);

%% fft3 for psf
FFT3PSFTranspose = fftForStack(rot90(largepsf,2));

% Delete variable 'largepsf' in the caller workspace as well as this workspace
% and create it again at the end of this function.
evalin('caller','clear largepsf'); clear largepsf;
% Although this may lead to tiny time waste, the released memory in GPU can
% be used in the next few steps, which gives fRL reconstruction ability for
% larger volmes.

%% backward projection for HXguess
HXguessBack = ifftForStack((repmat(HXguessFFT,[1,1,psf_s])) .* FFT3PSFTranspose);
% This line equals to
%     HXguessBackFFT = (repmat(HXguessFFT,[1,1,psf_s])) .* FFT3PSFTranspose;
%     HXguessBack = ifftForStack(HXguessBackFFT);
% but it seems it can reduce the maximum memory used.
        
%% backward projection for projNow
%%%%%%%%%%%%%%%
% This function is fully operated on GPU. Nevertheless, if an 
% out-of-gpu-memory error occurs here, you can gather 
% "projNowExpand" and "FFT3PSFTranspose" to move part of this backward
% projector on CPU.
% 
%     projNowExpand = zeros(proj_r,proj_c,psf_s,'single');
%     FFT3PSFTranspose = gather(FFT3PSFTranspose);
%%%%%%%%%%%%%%%

projNowExpand = gpuArray.zeros(proj_r,proj_c,psf_s,'single');
projNowExpand(:,:,(psf_s+1)/2) = projNow;

projBackfft = fftForStack(projNowExpand).* FFT3PSFTranspose;clear projNowExpand FFT3PSFTranspose;projBack=ifftForStack(projBackfft);clear projBackfft;
% This line equals to (both in grammer and in time-cost, but not in memory usage):
%     projBack = ifftForStack(fftForStack(projNowExpand).* FFT3PSFTranspose);%lines such as this is in great cost of memory
% This line also equals to:
%     projNowFFT = fftForStack(blurImageHereExpand);
%     projBackFFT = projNowFFT .* FFT3PSFTranspose;
%     projBack = ifftForStack(projBackFFT);

%% errorBack
errorBack = (real(projBack./HXguessBack)); % Calculate Error Matric

%% create largepsf again
evalin('caller','largepsf = gpuArray.zeros(proj_r,proj_c,psf_s,"single");');
end