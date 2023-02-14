%% Reconstruct captured projections
% PSFLoading - imageLoading - reconstruction
% ELi, 20220726
% ELi, 20230210, add multi-site DAO options

%% clear and paths
clc; clear; close all; addpath(genpath('utils'));
PSFs_Path = 'data//PSFs//25X'; %Check if PSF exists before generating new PSF
EI_path = 'data//images';%file path
EI_name = 'cx3cr1.tif';%file name
resultSave_path = ['results//recon_',datestr(now, 'YYYYmmDD_HHMMSS')];%Save results here

%% loading PSFs 
disp('Loading PSFs...');
for i = 1:13
    load([PSFs_Path,'//PSF_cx3cr1//psf_all_',num2str(i),'.mat'],'psf_thisAngle');
    psfs(:,:,:,i) = psf_thisAngle;
end
load([PSFs_Path,'//PSF_cx3cr1//PSFParameters.mat']);
disp('PSFs loaded');

%% loading images
disp('Loading raw data...');
projs = single(loadtiff([EI_path,'//',EI_name]));
disp('Raw data loaded');

%% reconstruct - DAO off
resultSave_path1 = [resultSave_path,'//DAO0']; % save results here
reconOpts.DAO = 0; % DAO off
%%%see reconPipeline.m for more reconstruct options
Xguess = reconPipeline(psfs,projs,PSFParameters,resultSave_path1,reconOpts);

%% reconstruct - DAO on
resultSave_path2 = [resultSave_path,'//DAO1']; % save results here
reconOpts.DAO = 1; % DAO on
[~,shiftMap] = reconPipeline(psfs,projs,PSFParameters,resultSave_path2,reconOpts);
phaseRecon(PSFParameters, shiftMap, resultSave_path2); % phase reconstruction

%% reconstruct - multi-site DAO
resultSave_path3 = [resultSave_path,'//DAO1_multiSite']; % save results here
reconOpts.DAO = 1; % DAO on
reconOpts.DAOOpts.patchN = 2; % (patchN * patchN) patches in total
[~,shiftMap] = reconPipeline(psfs,projs,PSFParameters,resultSave_path3,reconOpts);
phaseRecon(PSFParameters, shiftMap, resultSave_path3); % phase reconstruction