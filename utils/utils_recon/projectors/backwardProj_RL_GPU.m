% Backward projector in RL for 3D reconstruction
% 
% REFERENCE
%     Broxton M, Grosenick L, Yang S, et al. Wave optics theory and 3-D
%     deconvolution for the light field microscope [J]. Optics Express, 
%     2013, 21(21): 25418-25439.
% 
% ELi, 20230131, add comments
function projBack = backwardProj_RL_GPU(psf_t,proj)
% INPUT
%     psf_t   - transposed psf at some angle
%     proj    - projection at some angle

%% initialize&&size
[ra, ca]=size(proj);
[rb, cb, sb]=size(psf_t);
r = ra+rb;c=ca+cb; p1 = (r-ra)/2; p2=(c-ca)/2;
b1 = gpuArray.zeros(r,c,'single');
a1 = gpuArray.zeros(r,c,'single');
projBack=gpuArray.zeros(ra,ca,sb,'single');

%% convolve slice by slice
for z=1:size(psf_t,3)
    a1(1:ra,1:ca) = proj(:,:) ;
    b1(1:rb,1:cb) = psf_t(:,:,z) ;
    clear con1;
    con1 = ifft2(fft2(a1) .* fft2(b1));
    projBack(:,:,z) = real(con1(p1+1:r-p1,p2+1:c-p2));
end
end