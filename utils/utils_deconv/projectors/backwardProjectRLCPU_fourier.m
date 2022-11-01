function errorBack = backwardProjectRLCPU_fourier(largepsf,HXguessFFT,psf_s,img_r,img_c,blurImageHere)

        FFT3PSFTranspose = fftForStack(rot90(largepsf,2));
%         evalin('caller','clear largepsf');clear largepsf;%delete variable largepsf in the caller workspace and this workspace, create it again at the end of this function. Although this may lead to tiny time costs, the memory is saved.
        
        HXguessBack = ifftForStack((repmat(HXguessFFT,[1,1,psf_s])) .* FFT3PSFTranspose);
% This line equals to:
%         HXguessBackFFT = (repmat(HXguessFFT,[1,1,psf_s])) .* FFT3PSFTranspose;
%         HXguessBack = ifftForStack(HXguessBackFFT);
        

        %%%%%%%%%%%%%%%
        %if an out-of-gpu-memory occurs during this function, gather "blurImageHereExpand" "FFT3PSFTranspose" into cpu memory.
%         blurImageHereExpand = zeros(img_r,img_c,psf_s,'single');FFT3PSFTranspose = gather(FFT3PSFTranspose);
        %%%%%%%%%%%%%%%
        
        blurImageHereExpand = zeros(img_r,img_c,psf_s,'single');
        blurImageHereExpand(:,:,(psf_s+1)/2) = blurImageHere;

        acqBackfft = fftForStack(blurImageHereExpand).* FFT3PSFTranspose;clear blurImageHereExpand FFT3PSFTranspose;acqBack=ifftForStack(acqBackfft);clear acqBackfft;
% The above line equals to(both in grammer and in time-cost, but not in memory usage):
%         acqBack = ifftForStack(fftForStack(blurImageHereExpand).* FFT3PSFTranspose);%lines such as this is in great cost of memory
% This line also equals to:
%         acqBackFFT = acqFFT .* FFT3PSFTranspose;
%         acqBack = ifftForStack(acqBackFFT);

        errorBack = (real(acqBack./HXguessBack)); % Calculate Error Matric
%         evalin('caller','largepsf = gpuArray.zeros(img_r,img_c,psf_s,"single");');
end