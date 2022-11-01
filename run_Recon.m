%% Reconstruct captured images
% PSFLoading - imageLoading - reconstruction
% ELi, 20220726

%% clear and paths
clc; clear; close all; addpath(genpath('utils'));
PSFs_Path = 'data//PSFs//25X'; %Check if PSF exists before generating new PSF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EI_path = 'data//images';%file path
EI_name = 'cx3cr1.tif';%file name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultSave_path = ['results//recon_',datestr(now, 'YYYYmmDD_HHMMSS')];%Save results in this folder

%% loading PSFs 
disp('Loading PSFs...');
load([PSFs_Path,'//PSF_cx3cr1//psf_all.mat']);
load([PSFs_Path,'//PSF_cx3cr1//PSFParameters.mat']);
disp('PSFs loaded');

%% loading images
disp('Loading raw data...');
proj_all = single(loadtiff([EI_path,'//',EI_name]));
disp('Raw data loaded');
proj_all = preprocessing_stackResize(proj_all,[PSFParameters.EIParameters.upSample*size(proj_all,1),PSFParameters.EIParameters.upSample*size(proj_all,2)]); % upsample

%% reconstruct
deconvOptions.AO = 1; % DAO on
[~,shifts_estimated,deconvOptions] = deconvPipeline(psf_all,proj_all,PSFParameters,resultSave_path,deconvOptions);
disp('Reconstructing pupil phase...');
phaseReconstruction(PSFParameters,shifts_estimated,deconvOptions.deconvSavepath); % phase reconstruction
disp('Phase reconstructed.');
