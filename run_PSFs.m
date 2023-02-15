%% This script shows how to generate point spread functions (PSFs) in 2pSAM
% PSFGeneration
% ELi, 20220726

%% clear and paths
clc; clear; close all; addpath(genpath('utils'));
PSFs_savePath = ['results//PSFs_',datestr(now, 'YYYYmmDD_HHMMSS')];%Check if PSF exists before generating new PSF
systemParameters_path = 'data//systemParameters_simu//25X';

%% generate PSFs
%%%see PSFsGenerator.m for more PSF parameters
PSFParameters.imagingMode = 'mid';% use 'mid' for mid-NA 2pSAM configuration / use 'min' for min-NA 2pSAM configuration 
[psfs,PSFParameters] = PSFsGenerator(PSFParameters,systemParameters_path,PSFs_savePath);