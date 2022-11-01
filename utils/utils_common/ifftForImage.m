function output = ifftForImage(input)

if ndims(input)~=2
    error();
end

output = fftshift(ifft2(ifftshift(input)));