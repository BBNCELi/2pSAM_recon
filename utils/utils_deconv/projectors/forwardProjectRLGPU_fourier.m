function HXguessFFT = forwardProjectRLGPU_fourier(largepsf,Xguess,psf_s)

FFT3PSFFlip = fftForStack(flip(largepsf,3));
HXguessFFT = sum( fftForStack(Xguess) .* FFT3PSFFlip ,3)./psf_s;

end