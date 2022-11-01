% My self-coded ifftn function for an image stack.
% ELiiiiiii, 20201128

function output = ifftForStack(input)

if ndims(input)~=3
    if ndims(input) == 2
        warning('Using ifftForStack on an image...');
    else
        error('Dimention error when using ifftForStack!');
    end
end

output = fftshift(ifftn(ifftshift(input)));