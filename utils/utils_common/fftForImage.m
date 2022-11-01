function output = fftForImage(input)

if ndims(input)~=2
    error();
end

output = fftshift(fft2(ifftshift(input)));