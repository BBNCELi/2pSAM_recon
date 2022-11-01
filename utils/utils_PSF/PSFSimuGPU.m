% Prepare parameters and generate PSF in 2pSAM
% ELi, 20220530
function [psf,pupil,pupil_lowNA,pupil_ap] = PSFSimuGPU(aperture,shifty,shiftx,K,lambda, numAper, rindexObj, dxy, defocus, ...
    xysize, dkxy, kMaxPixels, lowNAkMaxPixels, pinholeDefocus, pupilPlaneMatrix, rindexSp, xysize_zoomFactor)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% preparations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gpuDevice;
psf = zeros(xysize_zoomFactor,xysize_zoomFactor,length(defocus),'single');
N = round((xysize+1)/2);

%Calculated the wavelength of light inside the objective lens and specimen
lambdaObj = lambda/rindexObj;
lambdaSp = lambda/rindexSp;
% Calculate the wave vectors in vaccuum, objective and specimens
k0 = 2*pi/lambda;
kObj = 2*pi/lambdaObj;
kSp = 2*pi/lambdaSp;
% Calculate pupil function
kxcord = (1:xysize)'-N;%Setting N as the center of the pupil
kycord = kxcord;
[kx, ky] = meshgrid(kxcord, kycord);
k = sqrt(kx.^2+ky.^2);
pupil = (k< kMaxPixels);
% Calculate pupil function in lowNA
kxcord = (1:xysize)'-N;%Setting N as the center of the pupil_lowNA
kycord = kxcord;
[kx, ky] = meshgrid(kxcord, kycord);
k = sqrt(kx.^2+ky.^2);
pupil_lowNA = (k< lowNAkMaxPixels);
% Calculate the sine of the semi-aperture angle in the objective lens
sinthetaObj = (k.*(dkxy))/kObj;
sinthetaObj(sinthetaObj>1) = 1;
% Calculate the cosine of the semi-aperture angle in the objective lens
costhetaObj = eps+sqrt(1-(sinthetaObj.^2));
% Calculate the sine of the semi-aperture angle in the specimen
sinthetaSp = (k.*(dkxy))/kSp;
sinthetaSp(sinthetaSp>1) = 1;
% Calculate the cosine of the semi-aperture angle in the specimen
costhetaSp = eps+sqrt(1-(sinthetaSp.^2));
% Defocus Phase calculation
phid = (sqrt(-1)*kObj).*costhetaObj;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spot just before pinhole
focal_p=ifftForImage(pupil_lowNA);
% spot just before pinhole, considering pinhole defocus
focal_pinhole = fresnelPropagationFor2DWavefront(focal_p, N, dkxy, pinholeDefocus, lambda, 1);
% spot just after pinhole, considering pinhole defocus
focal_afterPinhole=focal_pinhole.*aperture;
% spot just after pinhole
focal_ap = fresnelPropagationFor2DWavefront(focal_afterPinhole, N, dkxy, -pinholeDefocus, lambda, 1);
% K
focal_shiftK = fresnelPropagationFor2DWavefront(focal_ap, N, dkxy, K, lambda, 1);
% pupil, before angle scanning
pupil_p=fftForImage(focal_shiftK);
% pupil, after angle scanning
pupil_pshift=imtranslate(pupil_p,[shiftx,shifty]);
% pupil, considering the low-pass filtering of objective
pupil_ap=pupil_pshift.*pupil;
% pupil, considering aberrations
pupil_ap_withAberration = pupil_ap .* exp(-1i*pupilPlaneMatrix);
% defocus in z-axis
pupil_ap_withAberration = gpuArray(single(pupil_ap_withAberration) * length(pupil_ap_withAberration(:)));%To keep precision
zkBefore=0;
for zk = 1:length(defocus)
    zkCountingBlock = zk-zkBefore;
    OPDDefocus = exp(defocus(1, zk).*phid);
    pupilSA = pupil_ap_withAberration.*OPDDefocus;
    % Calculate coherent PSF using inverse Fourier Transform
    psf_temp(:, :, zkCountingBlock) = ifftForImage(pupilSA);
    if zkCountingBlock>4 || zk==length(defocus)
        psf_temp=gather(psf_temp);
        for imresizeCount = 1:zkCountingBlock
            % Calculate incoherent PSF from the coherent PSF
            % DownSample to save time
            psf(:,:,imresizeCount+zkBefore)=imresize(abs(psf_temp(:,:,imresizeCount)).^2,[xysize_zoomFactor,xysize_zoomFactor],'bilinear');
        end
        clear psf_temp;  
        zkBefore=zk;
    end
end