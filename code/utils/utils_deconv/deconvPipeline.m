% Deconvolution pipeline for 2pSAM
%%ELi, 20210721: Add auto-projections
%%ELi, 20220530: add shifts_estimated in output
function [Xguess,shifts_estimated,deconvOptions] = deconvPipeline(psf_all,proj_all,PSFParameters,resultSave_path,deconvOptions)
%% default % options
if ~exist('deconvOptions','var');deconvOptions = struct;end
if ~isfield(deconvOptions,'solver');deconvOptions.solver = 'fRL';end %'fRL_'/'fRL_avg'/'RL'/'fRL_avg_stdInit'
if ~isfield(deconvOptions,'maxIter');deconvOptions.maxIter = 10;end %0(for _simUpdate especially)/1/2/3/...
if ~isfield(deconvOptions,'maxIter_init');deconvOptions.maxIter_init = 1;end %1/2/3/...
if ~isfield(deconvOptions,'AO');deconvOptions.AO = 0;end %0/1/2
if ~isfield(deconvOptions,'equalWeight');deconvOptions.equalWeight = 0.5;end %0-1
if ~isfield(deconvOptions,'gpuSelect');deconvOptions.gpuSelect = 1;end %1/2/...
if ~isfield(deconvOptions,'startFrame');deconvOptions.startFrame = 1;end %1/14/27/...
if ~isfield(deconvOptions,'step');deconvOptions.step = PSFParameters.angleNum;end %13/26/39/...
if ~isfield(deconvOptions,'endFrame');deconvOptions.endFrame = size(proj_all,3);end %1300/100/...
deconvOptions.endFrame = min(deconvOptions.endFrame,size(proj_all,3));
if ~isfield(deconvOptions,'dropout');deconvOptions.dropout = 0;end %1/0/true/false
if ~isfield(deconvOptions,'itersTifSaveOn');deconvOptions.itersTifSaveOn = 0;end %1/0/true/false
if ~isfield(deconvOptions,'project_cutz');deconvOptions.project_cutz = 20;end %(bilateral)20/30/...
if ~isfield(deconvOptions,'project_cutx');deconvOptions.project_cutx = 20;end %(bilateral)20/30/...
if ~isfield(deconvOptions,'project_cuty');deconvOptions.project_cuty = 20;end %(bilateral)20/30/...
if ~isfield(deconvOptions,'demotion');deconvOptions.demotion = 0;end %0/1
%%%dir and save options
deconvOptions.deconvFolderName = deconvNameGenerator(deconvOptions,PSFParameters);
deconvOptions.deconvSavepath = [resultSave_path,'/',deconvOptions.deconvFolderName];
if ~exist(deconvOptions.deconvSavepath,'dir');mkdir(deconvOptions.deconvSavepath);end
save([deconvOptions.deconvSavepath,'/deconvOptions.mat'],'deconvOptions','PSFParameters','-v7.3');

%% preparations
psf_all = single(psf_all); proj_all = single(proj_all);
[~, ~, psfSize3ALSOsampleSize3, angleNum]=size(psf_all);
[imgSize1,imgSize2,imgSize3]=size(proj_all);
updateWeight=ones(angleNum,1)*deconvOptions.equalWeight;
% disp('Initializing GPU...');
% gpuDevice(deconvOptions.gpuSelect);
% disp('GPU initialization finished');
matSaveOn = 0;
tifSaveOn = deconvOptions.itersTifSaveOn;

%% Xguess initialization
switch deconvOptions.solver
    case {'fRL', 'fRL_avg', 'RL'}
        Xguess_init=ones(imgSize1,imgSize2,psfSize3ALSOsampleSize3);
        Xguess_init=Xguess_init./sum(Xguess_init(:))./angleNum;
    case {'fRL_avg_stdInit'}
        if isfield(deconvOptions,'Xguess_init')
            disp('Found Xguess_init...');
            Xguess_init = gpuArray(deconvOptions.Xguess_init);
        else
            disp('Calculating Xguess_init...');
            Xguess_init=ones(imgSize1,imgSize2,psfSize3ALSOsampleSize3);
            Xguess_init=Xguess_init./sum(Xguess_init(:))./angleNum;
            maxIter_stdInit = deconvOptions.maxIter;
            updateWeight_stdInit = updateWeight;
            if ~isfield(deconvOptions,'proj_all_init')
                proj_all_std = multiAngleStd(proj_all,angleNum);
                proj_all_std = proj_all_std - min(proj_all_std,[],[1,2]);
            else
                proj_all_std = deconvOptions.proj_all_init;
            end
            if deconvOptions.demotion
                proj_all_std = demotion(psf_all,proj_all_std,PSFParameters,deconvOptions);
            end
            Xguess_init = deconvRLGPU_fourier_DAO(psf_all,proj_all_std,PSFParameters,maxIter_stdInit,1,deconvOptions.AO,Xguess_init,updateWeight_stdInit,deconvOptions.deconvSavepath,matSaveOn,1);
        end
    otherwise
        error('Unexpected solver...!');
end

%% frame-by-frame deconvolution
for frameNow = deconvOptions.startFrame:deconvOptions.step:deconvOptions.endFrame
    if frameNow+angleNum-1 > imgSize3; break; end
%     disp(['frame ',num2str(frameNow),' reconstruction begins...']);
    disp(['Reconstruction begins...']);
    proj_all_forFrameNow = circshift(proj_all(:,:,frameNow:frameNow+angleNum-1),[0,0,frameNow-1]);
    psf_all_forFrameNow = psf_all;
    updateWeight_forFrameNow = updateWeight;
    %% demotion if needed
    if deconvOptions.demotion
        proj_all_forFrameNow = demotion(psf_all_forFrameNow,proj_all_forFrameNow,PSFParameters,deconvOptions);
    end
    %% select solver and strat deconvolution
    ticForThisFrame = tic;
    switch deconvOptions.solver
        case 'fRL'
            [Xguess,shifts_estimated] = deconvRLGPU_fourier_DAO(psf_all_forFrameNow,proj_all_forFrameNow,PSFParameters,deconvOptions.maxIter,deconvOptions.maxIter_init,deconvOptions.AO,Xguess_init,updateWeight_forFrameNow,deconvOptions.deconvSavepath,matSaveOn,tifSaveOn);
        case 'fRL_avg'
            [Xguess,shifts_estimated] = deconvRLGPU_fourier_DAO_avgUpdate(psf_all_forFrameNow,proj_all_forFrameNow,PSFParameters,deconvOptions.maxIter,deconvOptions.maxIter_init,deconvOptions.AO,Xguess_init,updateWeight_forFrameNow,deconvOptions.deconvSavepath,matSaveOn,tifSaveOn);
        case 'RL'
            [Xguess,shifts_estimated] = deconvRLGPU_DAO(psf_all_forFrameNow,proj_all_forFrameNow,PSFParameters,deconvOptions.maxIter,deconvOptions.maxIter_init,deconvOptions.AO,Xguess_init,updateWeight_forFrameNow,deconvOptions.deconvSavepath,matSaveOn,tifSaveOn);
        case 'fRL_avg_stdInit'
            [Xguess,shifts_estimated] = deconvRLGPU_fourier_DAO_avgUpdate(psf_all_forFrameNow,proj_all_forFrameNow,PSFParameters,0,1,0,Xguess_init,updateWeight_forFrameNow,deconvOptions.deconvSavepath,matSaveOn,tifSaveOn);
    end
%     disp(['frame ',num2str(frameNow),' reconstruction ends. Elapsed time: ',num2str(toc(ticForThisFrame)),' sec']);
%     disp(['Reconstruction ends. Elapsed time: ',num2str(toc(ticForThisFrame)),' sec']);
    disp(['Reconstruction ends']);

    %% volume save
    Xguess = gather(Xguess);
%     saveastiff_overwrite(Xguess, [deconvOptions.deconvSavepath,'/frame_',num2str(frameNow),'.tif']); 
    
    %% projs save - default
    thisFrame_xyzMiddle = Xguess(deconvOptions.project_cutx+1:end-deconvOptions.project_cutx,deconvOptions.project_cuty+1:end-deconvOptions.project_cuty,deconvOptions.project_cutz+1:end-deconvOptions.project_cutz);
    saveastiff_overwrite(thisFrame_xyzMiddle, [deconvOptions.deconvSavepath,'//output_3DStack',num2str(frameNow),'_dxy',num2str(PSFParameters.dxy/1000,'%.2f'),'dz',num2str(PSFParameters.defocusInterval),'um.tif']); 
    
    optionsForProjectionSave.append = true;
    optionsForProjectionSave.big = true;
    
%     gifFilename1 = [deconvOptions.deconvSavepath,'//MIPz-frame',num2str(deconvOptions.startFrame),'-.tif'];
%     thisFrame_xyzMiddle_mipz = max(thisFrame_xyzMiddle,[],3);
%     saveastiff(thisFrame_xyzMiddle_mipz,gifFilename1,optionsForProjectionSave);
%     gifFilename2 = [deconvOptions.deconvSavepath,'//MEANz-frame',num2str(deconvOptions.startFrame),'-.tif'];
%     thisFrame_xyzMiddle_meanz = mean(thisFrame_xyzMiddle,3);
%     saveastiff(thisFrame_xyzMiddle_meanz,gifFilename2,optionsForProjectionSave);
%     gifFilename4 = [deconvOptions.deconvSavepath,'//MIPy-frame',num2str(deconvOptions.startFrame),'-.tif'];
%     thisFrame_xyzMiddle_mipy = squeeze(max(thisFrame_xyzMiddle,[],2));
%     saveastiff(thisFrame_xyzMiddle_mipy,gifFilename4,optionsForProjectionSave);
%     gifFilename5 = [deconvOptions.deconvSavepath,'//MEANy-frame',num2str(deconvOptions.startFrame),'-.tif'];
%     thisFrame_xyzMiddle_meany = squeeze(mean(thisFrame_xyzMiddle,2));
%     saveastiff(thisFrame_xyzMiddle_meany,gifFilename5,optionsForProjectionSave);
%     gifFilename6 = [deconvOptions.deconvSavepath,'//MIPx-frame',num2str(deconvOptions.startFrame),'-.tif'];
%     thisFrame_xyzMiddle_mipx = squeeze(max(thisFrame_xyzMiddle,[],1));
%     saveastiff(thisFrame_xyzMiddle_mipx,gifFilename6,optionsForProjectionSave);
%     gifFilename7 = [deconvOptions.deconvSavepath,'//MEANx-frame',num2str(deconvOptions.startFrame),'-.tif'];
%     thisFrame_xyzMiddle_meanx = squeeze(mean(thisFrame_xyzMiddle,1));
%     saveastiff(thisFrame_xyzMiddle_meanx,gifFilename7,optionsForProjectionSave);


    gifFilename1 = [deconvOptions.deconvSavepath,'//output_3DStack_z.tif'];
    thisFrame_xyzMiddle_mipz = max(thisFrame_xyzMiddle,[],3);
    thisFrame_xyzMiddle_mipz = uint8(thisFrame_xyzMiddle_mipz / (0.3*max(thisFrame_xyzMiddle_mipz(:)))*255);
%     saveastiff(thisFrame_xyzMiddle_mipz,gifFilename1,optionsForProjectionSave);
        imwrite(thisFrame_xyzMiddle_mipz,gifFilename1);
    gifFilename4 = [deconvOptions.deconvSavepath,'//output_3DStack_y.tif'];
    thisFrame_xyzMiddle_mipy = squeeze(max(thisFrame_xyzMiddle,[],2));
    thisFrame_xyzMiddle_mipy = uint8(thisFrame_xyzMiddle_mipy / (0.7*max(thisFrame_xyzMiddle_mipy(:)))*255);
%     saveastiff(thisFrame_xyzMiddle_mipy,gifFilename4,optionsForProjectionSave);
        imwrite(thisFrame_xyzMiddle_mipy,gifFilename4);
    gifFilename6 = [deconvOptions.deconvSavepath,'//output_3DStack_x.tif'];
    thisFrame_xyzMiddle_mipx = squeeze(max(thisFrame_xyzMiddle,[],1));
    thisFrame_xyzMiddle_mipx = uint8(thisFrame_xyzMiddle_mipx / (0.7*max(thisFrame_xyzMiddle_mipx(:)))*255);
%     saveastiff(thisFrame_xyzMiddle_mipx,gifFilename6,optionsForProjectionSave);
        imwrite(thisFrame_xyzMiddle_mipx,gifFilename6);


    %% projs save - optional
    if isfield(deconvOptions,'project_opt')
        for project_opt_count = 1:length(deconvOptions.project_opt)
            project_opt_here = deconvOptions.project_opt{project_opt_count};
            project_opt_here_fileName = [deconvOptions.deconvSavepath,'/',project_opt_here{1},'_',num2str(project_opt_here{2}),'_',num2str(project_opt_here{3}),'_',num2str(project_opt_here{4}),'_',num2str(project_opt_here{5}),'_',num2str(project_opt_here{6}),'_',num2str(project_opt_here{7}),'-frame',num2str(deconvOptions.startFrame),'-.tif'];
            switch project_opt_here{1}(1:end-1)
                case 'MIP'
                    switch project_opt_here{1}(end)
                        case 'x'
                            project_result_here = squeeze(max(Xguess(project_opt_here{2}:project_opt_here{3},project_opt_here{4}:project_opt_here{5},project_opt_here{6}:project_opt_here{7}),[],1));
                        case 'y'
                            project_result_here = squeeze(max(Xguess(project_opt_here{2}:project_opt_here{3},project_opt_here{4}:project_opt_here{5},project_opt_here{6}:project_opt_here{7}),[],2));
                        case 'z'
                            project_result_here = squeeze(max(Xguess(project_opt_here{2}:project_opt_here{3},project_opt_here{4}:project_opt_here{5},project_opt_here{6}:project_opt_here{7}),[],3));
                    end
                case 'MEAN'
                    switch project_opt_here{1}(end)
                        case 'x'
                            project_result_here = squeeze(mean(Xguess(project_opt_here{2}:project_opt_here{3},project_opt_here{4}:project_opt_here{5},project_opt_here{6}:project_opt_here{7}),1));
                        case 'y'
                            project_result_here = squeeze(mean(Xguess(project_opt_here{2}:project_opt_here{3},project_opt_here{4}:project_opt_here{5},project_opt_here{6}:project_opt_here{7}),2));
                        case 'z'
                            project_result_here = squeeze(mean(Xguess(project_opt_here{2}:project_opt_here{3},project_opt_here{4}:project_opt_here{5},project_opt_here{6}:project_opt_here{7}),3));
                    end
            end
            saveastiff(project_result_here,project_opt_here_fileName,optionsForProjectionSave);
        end
    end
end
