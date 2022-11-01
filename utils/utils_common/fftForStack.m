% My self-coded fftn function for an image stack.
% ELiiiiiii, 20201128

function output = fftForStack(input)

if ndims(input)~=3
    if ndims(input) == 2
        warning('Using fftForStack on an image...');
    else
        error('Dimention error when using fftForStack!');
    end
end

output = fftshift(fftn(ifftshift(input)));