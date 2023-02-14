% Forward projector in RL for 3D reconstruction
% 
% REFERENCE
%     Broxton M, Grosenick L, Yang S, et al. Wave optics theory and 3-D
%     deconvolution for the light field microscope [J]. Optics Express, 
%     2013, 21(21): 25418-25439.
% 
% ELi, 20230131, add comments
function projection = forwardProj_RL_GPU(psf,Xguess)
% INPUT
%     psf       - psf at some angle
%     Xguess    - volume given

%% initialize&&size
projection=gpuArray.zeros(size(Xguess,1),size(Xguess,2),'single');
[ra, ca, sa]=size(Xguess);
[rb, cb, sb]=size(psf);
r = ra+rb; c=ca+cb; p1 = (r-ra)/2; p2=(c-ca)/2;
a1 = gpuArray.zeros(r,c,'single');
b1 = gpuArray.zeros(r,c,'single');

%% convolve slice by slice and sum
for z=1:size(Xguess,3)
    a1(1:ra,1:ca) = Xguess(:,:,z) ;
    b1(1:rb,1:cb) = psf(:,:,z) ;
    clear con1;
    con1 = ifft2(fft2(a1) .* fft2(b1));
    projection = projection + real(con1(p1+1:r-p1,p2+1:c-p2));
end
end