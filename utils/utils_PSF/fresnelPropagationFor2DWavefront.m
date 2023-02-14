% The old fresnel2D function is a mess.
% Rewrite this function strictly according to Goodman's book.

% Inputs:
%         inputWavefront: the WAVEFRONT of light at z = 0nm. Should include amplitude AND phase.
%         centrePixelIndex: setting inputWavefront(centrePixelIndex,centrePixelIndex) as the centre of inputWavefront.
%         dkxy: The frequence distance between pixels. 
%             USUALLY we use dkxy = (2*pi)/(xysize*dxy) 
%                 where xysize is the number of columns (== the number of rows) of inputWavefront
%                       dxy is the physical distance between pixels.
%         z: nm. The expected distance between outputWavefront and inputWavefront.
%         lambda: nm. The wavelength in vacuum.
%         rIndex: Refractive index of the medium between outputWavefront and inputWavefront.
% Outputs:
%         outputWavefront: the WAVEFRONT of light at z = 'z_nm' nm. Should include amplitude AND phase.

% ELiiiiiii, 20201117
function outputWavefront=fresnelPropagationFor2DWavefront(inputWavefront,centrePixelIndex,dkxy,z_nm,lambda_nm,rIndex)

[xsize,ysize] = size(inputWavefront);
if xsize~= ysize
    error('The inputWavefront of function "fresnelPropagationFor2DWavefront" should be a square!!!');
end

%Calculated the wavelength of light inside this medium
lambdaInThisMedium = lambda_nm/rIndex;

% Calculate the wave vectors in this medium, 
kInThisMedium = 2*pi/lambdaInThisMedium;

% Generate the pupil function amplitude
kxcord = ((1:xsize)'-centrePixelIndex)*dkxy; %Setting centrePixelIndex as the center of the pupil
kycord = kxcord;
[kx, ky] = meshgrid(kxcord, kycord);
k = sqrt(kx.^2+ky.^2);

% normalized and calculate kz
sqrtkxSquarePluskySquare = k/kInThisMedium;
sqrtkxSquarePluskySquare(sqrtkxSquarePluskySquare>1) = 1;
% kz = sqrt(1-(sqrtkxSquarePluskySquare.^2));
kz = 1-1/2*sqrtkxSquarePluskySquare.^2;%This is actually the fresnel diffraction

% Defocus Phase calculation
phid = (1i*kInThisMedium).*kz;
OPDDefocus = exp(z_nm.*phid);

% Calculate output wavefront
outputWavefront = ifftForImage(fftForImage(inputWavefront).*OPDDefocus);

