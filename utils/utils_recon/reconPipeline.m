% Reconstruction pipeline for 2pSAM
%
% ELi, 20230208, rearrange inputs, delete unused functions and add comments
function [Xguess,shiftMap,dispMap,reconOpts] = reconPipeline(psfs,projs,PSFParameters,saveFolder,reconOpts)
% INPUT
%     psfs      - 3D psfs used for reconstruction, (usually) 4D matrix
%     projs     - 2D projections to be reconstructed, (usually) 3D matrix
%     PSFParameters - the parameters of PSFs used for imaging
%     saveFolder - save reconstruction results&&options in this folder
%     reconOpts - options for 3D reconstruction, struct (default)
%                 .maxIter: max iterations (10)
%                 .solver: reconstruction algorithm (fRL)
%                 .upWeight: weight for volume update in each iteration (.5)
%                 .upSeq: updating sequence for projections from different angles
%                 .CAIndex: index of centre projection in projs (1)
%                 .initMethod: volume initialization method ('all1')
%                 .gpuSelect: GPU used (1)
%                 .startFrame: frame to start reconstruction (1)
%                 .step: skip frames to reconstruction (13)
%                 .endFrame: frame to end reconstruction (inf)
%                 .savepath: folder to save results. If not given, generate a folder in saveFolder and save results there
%                 .tifSaveOn: iterations save flag (.tif), 0/1 (0)
%                 .matSaveOn: iterations save flag (.mat), 0/1 (0)
%                 .dispOn: display reconstruction process flag, 0/1 (1)
%                 .project_cutz: slices cut for projection (bilateral) (20)
%                 .project_cutx: pixels cut in x-dim for projection (bilateral) (0)
%                 .project_cuty: pixels cut in y-dim for projection (bilateral) (0)
%                 .project_opt: self defined project options (e.g. {'MEANy',125,130,506,112,522,151})
%                 .DAO: DAO flag, 0/1
%                 .DAOOpts: options for DAO
%                     .maxIter_init: max iteration without DAO
%                     .shiftEstFunc: function used for shift estimation
%                     .patchN: number of patches for multi-site shift estimation
%                     .patchOvFactor: overlap factor for multi-site shift estimation
%                     .minPatchSize: minimum patch size for multi-site shift estimation
%                     .sidelobe: sidelobe in pixel for shift estimation
%                     .sidelobez: use only the high-resolution regime for shift estimation
%                 .demotion: demotion flag, 0/1
%                 .demotionOpts: options for motion correction
%                     .maxIter_init: iterations without shift estimation in each iteration of demotion
%                     .maxIter: iterations without+out shift estimation in each iteration of demotion
%                     .shiftEstFunc: function used for shift estimation
%                     .patchN: number of patches for multi-site shift estimation
%                     .patchOvFactor: overlap factor for multi-site shift estimation
%                     .minPatchSize: minimum patch size for multi-site shift estimation
%                     .sidelobe: sidelobe in pixel for shift estimation
%                     .sidelobez: use only the high-resolution regime for shift estimation
%                     .thrMSE: the motion map converges if its increments are small enough
%                     .maxIterDemotion: max iterations to quit demotion if the motion map does not converge
%                     .maxShift: max shift in pixel for each iteration of demotion
%                     .angDistx_normed: PSF angle distribution, x, to remove the defocus item when initializing motion map
%                     .angDisty_normed: PSF angle distribution, y, to remove the defocus item when initializing motion map

%% default recon options
if ~exist('reconOpts','var');reconOpts = struct;end
% .maxIter: max iterations (10)
if ~isfield(reconOpts,'maxIter');reconOpts.maxIter = 10;end
% .solver: reconstruction algorithm (fRL): fRL/RL
if ~isfield(reconOpts,'solver');reconOpts.solver = 'fRL';end
% .upWeight: weight for volume update in each iteration (.5): 0-1
if ~isfield(reconOpts,'upWeight');reconOpts.upWeight = 0.5;end
% .upSeq: updating sequence for projections from different angles
% .CAIndex: index of centre projection in projs (1)
[upSeq_def,CAInd_def] = defaultUpSeq(PSFParameters);
if ~isfield(reconOpts,'upSeq');reconOpts.upSeq = upSeq_def;end
if ~isfield(reconOpts,'CAIndex');reconOpts.CAIndex = CAInd_def;end
% .initMethod: volume initialization method ('all1'): 'all1'/'std'
% % 'std': use time-std to initialize volume when long-time projections need to be reconstructed
if ~isfield(reconOpts,'initMethod');reconOpts.initMethod = 'all1';end
% .gpuSelect: GPU used (1)
if ~isfield(reconOpts,'gpuSelect');reconOpts.gpuSelect = 1;end
% .startFrame: frame to start reconstruction (1)
if ~isfield(reconOpts,'startFrame');reconOpts.startFrame = 1;end
% .step: skip frames to reconstruction (13)
if ~isfield(reconOpts,'step');reconOpts.step = PSFParameters.angleNum;end
% .endFrame: frame to end reconstruction (inf)
if ~isfield(reconOpts,'endFrame');reconOpts.endFrame = size(projs,3);end
reconOpts.endFrame = min(reconOpts.endFrame,size(projs,3));
% .tifSaveOn: iterations save flag (.tif), 0/1 (0)
if ~isfield(reconOpts,'tifSaveOn');reconOpts.tifSaveOn = 1;end
% .matSaveOn: iterations save flag (.mat), 0/1 (0)
if ~isfield(reconOpts,'matSaveOn');reconOpts.matSaveOn = 0;end
% .dispOn: display reconstruction process flag, 0/1 (1)
if ~isfield(reconOpts,'dispOn');reconOpts.dispOn = 1;end
% .project_cutz: slices cut for projection (bilateral) (20)
if ~isfield(reconOpts,'project_cutz');reconOpts.project_cutz = 20;end %(bilateral)20/30/...
% .project_cutx: pixels cut in x-dim for projection (bilateral) (0)
if ~isfield(reconOpts,'project_cutx');reconOpts.project_cutx = 0;end %(bilateral)20/30/...
% .project_cuty: pixels cut in y-dim for projection (bilateral) (0)
if ~isfield(reconOpts,'project_cuty');reconOpts.project_cuty = 0;end %(bilateral)20/30/...
% .DAO: DAO flag, 0/1
if ~isfield(reconOpts,'DAO');reconOpts.DAO = 0;end
% .demotion: demotion flag, 0/1
if ~isfield(reconOpts,'demotion');reconOpts.demotion = 0;end
if reconOpts.demotion % the demotion algorithm need to know the angle distribution in the back pupil plane
    reconOpts.demotionOpts.angDistx_normed = PSFParameters.angDistx_normed;
    reconOpts.demotionOpts.angDisty_normed = PSFParameters.angDisty_normed;
end
% .savepath: folder to save results. If not given, generate a folder in saveFolder and save results there
savepath_def = [saveFolder,'//',saveFolderGenerator(reconOpts,PSFParameters)];
if ~isfield(reconOpts,'savepath');reconOpts.savepath = savepath_def;end
if ~exist(reconOpts.savepath,'dir');mkdir(reconOpts.savepath);end
%%% save recon options in savepath
save([reconOpts.savepath,'//reconOpts.mat'],'reconOpts','PSFParameters','-v7.3');

%% size and other preparations
psfs = single(psfs); projs = single(projs); % singular resolution is enough in most cases
[~, ~, ~, angleNum]=size(psfs); [~,~,frameNum]=size(projs);
gpuDevice(reconOpts.gpuSelect);

%% Xguess initialization
[Xguess_init, reconOpts] = initXguess(psfs,projs,reconOpts.initMethod,reconOpts);

%% frame-by-frame reconstruction
for frameNow = reconOpts.startFrame:reconOpts.step:reconOpts.endFrame
    if frameNow+angleNum-1 > frameNum; break; end

    disp(['frame ',num2str(frameNow),' reconstruction begins...']);
    projs_frameNow = projs(:,:,frameNow:frameNow+angleNum-1);
    psfs_frameNow = circshift(psfs,[0,0,0,1-frameNow]);

    %% demotion if needed
    if reconOpts.demotion
        projs_frameNow = demotion(psfs_frameNow,projs_frameNow,reconOpts);
    end

    %% select solver and start reconstruction
    ticForThisFrame = tic;

    switch reconOpts.solver
        case 'fRL'
            [Xguess,shiftMap,dispMap] = recon_fRL_GPU(psfs_frameNow,projs_frameNow,reconOpts,Xguess_init);
        case 'RL'
            [Xguess,shiftMap,dispMap] = recon_RL_GPU(psfs_frameNow,projs_frameNow,reconOpts,Xguess_init);
    end

    disp(['frame ',num2str(frameNow),' reconstruction ends. Elapsed time: ',num2str(toc(ticForThisFrame)),' sec']);

    %% save volume and shiftMap
    Xguess = gather(Xguess);
    saveastiff_overwrite(Xguess, [reconOpts.savepath,'/frame_',num2str(frameNow),'.tif']);
    %%%turn it on if needed
    save([reconOpts.savepath,'/frame_',num2str(frameNow),'_DAO.mat'], 'shiftMap', 'dispMap');
    
    %% project volume and save a video - default
    Xguess_xyzMiddle = Xguess(reconOpts.project_cutx+1:end-reconOpts.project_cutx,reconOpts.project_cuty+1:end-reconOpts.project_cuty,reconOpts.project_cutz+1:end-reconOpts.project_cutz);
    saveOpts.append = true;
    saveOpts.big = true;
    
    %%%MIPz
    gifFilename1 = [reconOpts.savepath,'/MIPz-frame',num2str(reconOpts.startFrame),'-.tif'];
    thisFrame_xyzMiddle_mipz = max(Xguess_xyzMiddle,[],3);
    saveastiff(thisFrame_xyzMiddle_mipz,gifFilename1,saveOpts);
    %%%MEANz
    gifFilename2 = [reconOpts.savepath,'/MEANz-frame',num2str(reconOpts.startFrame),'-.tif'];
    thisFrame_xyzMiddle_meanz = mean(Xguess_xyzMiddle,3);
    saveastiff(thisFrame_xyzMiddle_meanz,gifFilename2,saveOpts);
    %%%color depth, turn it on if needed
%    gifFilename3 = [deconvOptions.deconvSavepath,'/colorDepthz-frame',num2str(deconvOptions.startFrame),'-.tif'];
%    saveastiff_colorDepth(thisFrame_xyzMiddle,gifFilename3,optionsForProjectionSave);
    %%%MIPy
    gifFilename4 = [reconOpts.savepath,'/MIPy-frame',num2str(reconOpts.startFrame),'-.tif'];
    thisFrame_xyzMiddle_mipy = squeeze(max(Xguess_xyzMiddle,[],2));
    saveastiff(thisFrame_xyzMiddle_mipy,gifFilename4,saveOpts);
    %%%MEANy
    gifFilename5 = [reconOpts.savepath,'/MEANy-frame',num2str(reconOpts.startFrame),'-.tif'];
    thisFrame_xyzMiddle_meany = squeeze(mean(Xguess_xyzMiddle,2));
    saveastiff(thisFrame_xyzMiddle_meany,gifFilename5,saveOpts);
    %%%MIPx
    gifFilename6 = [reconOpts.savepath,'/MIPx-frame',num2str(reconOpts.startFrame),'-.tif'];
    thisFrame_xyzMiddle_mipx = squeeze(max(Xguess_xyzMiddle,[],1));
    saveastiff(thisFrame_xyzMiddle_mipx,gifFilename6,saveOpts);
    %%%MEANx
    gifFilename7 = [reconOpts.savepath,'/MEANx-frame',num2str(reconOpts.startFrame),'-.tif'];
    thisFrame_xyzMiddle_meanx = squeeze(mean(Xguess_xyzMiddle,1));
    saveastiff(thisFrame_xyzMiddle_meanx,gifFilename7,saveOpts);
    
    %% project volume and save a video - optional
    if isfield(reconOpts,'project_opt')
        for project_opt_count = 1:length(reconOpts.project_opt)
            project_opt_here = reconOpts.project_opt{project_opt_count};
            project_opt_here_fileName = [reconOpts.deconvSavepath,'/',project_opt_here{1},'_',num2str(project_opt_here{2}),'_',num2str(project_opt_here{3}),'_',num2str(project_opt_here{4}),'_',num2str(project_opt_here{5}),'_',num2str(project_opt_here{6}),'_',num2str(project_opt_here{7}),'-frame',num2str(reconOpts.startFrame),'-.tif'];
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
            saveastiff(project_result_here,project_opt_here_fileName,saveOpts);
        end
    end
end
