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

%% reconstruct
%%%see reconPipeline.m for more options
reconOpts.DAO = 1; % DAO on
[~,shiftMap] = reconPipeline(psfs,projs,PSFParameters,resultSave_path,reconOpts);
phaseRecon(PSFParameters, shiftMap, resultSave_path); % phase reconstruction