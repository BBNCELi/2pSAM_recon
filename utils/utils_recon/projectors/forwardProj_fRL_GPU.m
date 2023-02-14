% Forward projector in fRL for 3D reconstruction
% 
% REFERENCE
%     Wu J, et al. Iterative tomography with digital adaptive optics permits
%     hour-long intravital observation of 3D subcellular dynamics at millisecond
%     scale [J]. Cell, 2021, 184(12): 3318-3332 e3317.
% 
% ELi, 20230131, add comments
function HXguessFFT = forwardProj_fRL_GPU(largepsf,Xguess)
% INPUT
%     largepsf  - psf but in the size of projection
%     Xguess    - volume given

%% size
psf_s = size(largepsf,3);

%% fft3 for psf
FFT3PSFFlip = fftForStack(flip(largepsf,3));

%% forward projection
HXguessFFT = sum( fftForStack(Xguess) .* FFT3PSFFlip ,3)./psf_s;

end