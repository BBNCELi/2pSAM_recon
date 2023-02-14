% Initialize Xguess before 3D reconstruction
% 
% NOTES
% 
% ELi, 20230208
function [Xguess_init,reconOpts] = initXguess(psfs,projs,method,reconOpts)
% INPUT
%     psfs    - 3D psfs used for reconstruction, (usually) 4D matrix
%     projs   - 2D projections to be reconstructed, (usually) 3D matrix
%     method  - initialization method:
%               'all1': all 1 matrix (default)
%               'std': calculate the time-std of projs and reconstruct a volume as initialization

%% check inputs
%%%size and default
[psf_r, psf_c, psf_s, angleNum]=size(psfs);
[proj_r,proj_c,proj_num]=size(projs);
if ~exist('method','var')||isempty(method); method = 'all1'; end

switch method
%% case 'all1'
    case 'all1'
        Xguess_init = ones(proj_r,proj_c,psf_s);
        Xguess_init = Xguess_init./sum(Xguess_init(:))./angleNum;

%% case 'std'
    case 'std'
        disp('Initiating volume with time-std of projections...');
        projs_std = multiAngleStd(projs,angleNum);
        projs_std = projs_std - min(projs_std,[],[1,2]);
        if reconOpts.demotion
            projs_std = demotion(psfs,projs_std,reconOpts);

        end
        Xguess_init = ones(proj_r,proj_c,psf_s);
        Xguess_init = Xguess_init./sum(Xguess_init(:))./angleNum;
        Xguess_init = recon_fRL_GPU(psfs,projs_std,reconOpts,Xguess_init);
        % As std-reconstructed volume has incorporated the spatial information,
        % further reconstructions need only 1 iteration.
        reconOpts.maxIter = 1;
        reconOpts.demotion = 0;

%% otherwise
    otherwise
        error('Unexpected initialization method');
end