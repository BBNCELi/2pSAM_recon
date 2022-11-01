%Preprocessing pipeline for 2pSAM
%Parameters set ---> preprocessing ---> save results

%%ELi, 20210719
%%ELi, 20210725
%%ELi, 20210916: add 3D median filter to medianFilterSize
%%ELi, 20211101: add group-mean; add meanFilter3; add gaussianFilter3
%%ELi, 20220531: final version
function [proj_all, prepParameters, returnParameters] = preprocessPipeline(proj_all, angleNum, resultSave_path, varargin)
% varargin:
%%%%%% name %%%%%%%%%%%%%%% meaning %%%%%%%%% input examples %%%%%
%     'save'             -- on/off         -- 1/0/true/false
%     'dropout'          -- on/off         -- 1/0/true/false
%     'motionCorr'       -- motionCorrStar -- 1/true/0/false/512*512*13matrix
%     'removeOutOfPupil' -- zoom factor    -- 1.3/1.5/3/8/...
%     'deepcad'          -- model          -- 'model1'/"model1"/{'model1','model2','model3'}/...
%     'medianFilterSize' -- patch size     -- 3/7/[3,3,1]/[9,9,3]/[7,3]...
%     'meanFilterSize'   -- patch size     -- 3/7/[3,3,1]/[9,9,3]/[7,3]...
%     'gaussFilterSize'  -- sigma value    -- 3/7/[3,3,1]/[9,9,3]/[7,3]...
%     'backSubtract'     -- options        -- optionsForBackSubtract
%     'upsample2'        -- pixels         -- 512/1024/2048/...
%     'groupMean'        -- groupSize      -- 2/5/10/...
%% Parameters set
p = inputParser;
%%% save or not
p.addOptional(... 
  'save', 0, ... 
   @(x) validateattributes(x, {'numeric', 'logical'},...
   {'scalar'}));
%%% dropout or not
p.addOptional(...
    'dropout',1,...
    @(x) validateattributes(x, {'numeric', 'logical'},...
   {'scalar'}));
%%% motion correction or not --> motionCorrStar
p.addOptional(... 
  'motionCorr', 0, ... 
   @(x) validateattributes(x, {'numeric', 'logical'},...
   {}));
%%% deepcad or not --> model name(s)
p.addOptional(...
    'deepcad',-1,...
    @(x) validateattributes(x, {'string','char','cell'},...
    {}));
%%% remove out of pupil or not --> zoom factor
p.addOptional(...
    'removeOutOfPupil',-1,...
    @(x) validateattributes(x, {'double'},...
    {'scalar','positive'}));
%%% median filter or not --> patch size of median filter
p.addOptional(...
    'medianFilterSize',-1,...
    @(x) validateattributes(x, {'numeric'},...
    {'vector', 'integer', 'nonnegative'}));
%%% mean filter or not --> patch size of mean filter
p.addOptional(...
    'meanFilterSize',-1,...
    @(x) validateattributes(x, {'numeric'},...
    {'vector', 'integer', 'nonnegative'}));
%%% gaussian filter or not --> sigma value of gaussian filter
p.addOptional(...
    'gaussFilterSize',-1,...
    @(x) validateattributes(x, {'numeric'},...
    {'vector', 'integer', 'nonnegative'}));
%%% backSubtract or not --> optionsOfBackSubtract
p.addOptional(...
    'backSubtract',-1,...
    @(x) isa(x, 'struct'));
%%% upsampling or not --> upsample to
p.addOptional(...
    'upsample2',-1,...
    @(x) validateattributes(x, {'numeric'},...
    {'scalar', 'integer', 'positive'}));
%%% groupMean or not --> groupSize
p.addOptional(...
    'groupMean',-1,...
    @(x) validateattributes(x, {'numeric'},...
    {'scalar', 'integer', 'positive'}));
parse(p,varargin{:});

%% preparations
%%% save folders
if ~exist(resultSave_path,'dir');mkdir(resultSave_path);end

%% load existed results if there is one
potentialResultList = findNameList(resultSave_path,'preprocess','_parser.mat');
existedFileFound = 0;
for subFileCount = 1:length(potentialResultList)
    potentialResultName = [resultSave_path,'/',potentialResultList{subFileCount}];
    load(potentialResultName,'prepParameters');%load existed prepParameters
    if isequal(p.Results, prepParameters)
        disp('Found existed preprocessing results... Loading now...');
        load([potentialResultName(1:end-11),'.mat'],'proj_all');%load existed conv_sum
        load(potentialResultName,'returnParameters');%load existed returnParameters
        existedFileFound = 1;
        break;
    end
end

%% if not, process data according to given parameters
timeNow = datestr(now, 'YYYYmmDD_HHMMSS');
optionsForTiffSave.overwrite = true;
if existedFileFound ~= 1
    returnParameters = struct;
    saveOn = p.Results.save;
    %%% sub save folders
    if saveOn ~= 0
        folderNamePre_preprocessing = [resultSave_path,'/prep_',timeNow,'_angle'];
        folderName_preprocessing = cell(angleNum,1);
        for angleCount = 1:angleNum
            folderName_preprocessing{angleCount} = [folderNamePre_preprocessing,num2str(angleCount)];
            mkdir(folderName_preprocessing{angleCount});
        end
    end

    %%% dropout
    if p.Results.dropout == 1
        dropOutFrames = mod(size(proj_all,3), angleNum);
        proj_all = proj_all(:,:,1:end-dropOutFrames);
    end

    %%% motion correction
    returnParameters.motionCorr = 0;
    if ~isequal(p.Results.motionCorr, 0)
        expectedTemplateSize = [size(proj_all,1),size(proj_all,2),angleNum];
        if ~isequal(p.Results.motionCorr, 1) && ~isequal(size(p.Results.motionCorr),expectedTemplateSize)
            error('Unrecognized motionCorr input!');
        end
        returnParameters.motionCorr = zeros(expectedTemplateSize,'like',proj_all);
        for angleCount = 1:angleNum
            timeSequence_thisAngle = proj_all(:,:,angleCount:angleNum:end);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/beforeMotionCorr.tiff'],optionsForTiffSave);
            end
            if isequal(p.Results.motionCorr, 1)%%estimate templates automatically
                [timeSequence_thisAngle, templateHere] = preprocessing_2DRegistration(timeSequence_thisAngle);
            elseif isequal(size(p.Results.motionCorr),expectedTemplateSize)%%use input templates
                [timeSequence_thisAngle, templateHere] = preprocessing_2DRegistration(timeSequence_thisAngle,p.Results.motionCorr(:,:,angleCount));
            end
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/afterMotionCorr.tiff'],optionsForTiffSave);
            end
            proj_all(:,:,angleCount:angleNum:end) = timeSequence_thisAngle;
            returnParameters.motionCorr(:,:,angleCount) = templateHere;
        end
    end

    %%% deepcad
    if ~isequal(p.Results.deepcad, -1)
        models = p.Results.deepcad;
        if isa(models, 'string') || isa(models, 'char')
            models = repelem({models},angleNum);
        elseif isa(models, 'cell')
            if length(models) ~= angleNum
                error('Mismatch between input models and angle number.');
            end
        end
        for angleCount = 1:angleNum
            timeSequence_thisAngle = proj_all(:,:,angleCount:angleNum:end);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/beforeDeepcad.tiff'],optionsForTiffSave);
            end
            timeSequence_thisAngle = preprocessing_denoise_deepcad(timeSequence_thisAngle,char(models{angleCount}));
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/afterDeepcad.tiff'],optionsForTiffSave);
            end
            proj_all(:,:,angleCount:angleNum:end) = timeSequence_thisAngle;
        end
%         disp('Deepcad in Matlab currently not supported...');
    end

    %%% removeOutOfPupil
    if ~isequal(p.Results.removeOutOfPupil, -1) && ~isequal(p.Results.removeOutOfPupil, 0)
        [proj_all,~] = removeOutOfPupilHighFrequency(proj_all,p.Results.removeOutOfPupil);
    end

    %%% median filter
    if ~isequal(p.Results.medianFilterSize, -1)
        for angleCount = 1:angleNum
            timeSequence_thisAngle = proj_all(:,:,angleCount:angleNum:end);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/beforeMedianFilter.tiff'],optionsForTiffSave);
            end
            timeSequence_thisAngle = preprocessing_denoise_medianFilter3(timeSequence_thisAngle,p.Results.medianFilterSize);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/afterMedianFilter.tiff'],optionsForTiffSave);
            end
            proj_all(:,:,angleCount:angleNum:end) = timeSequence_thisAngle;
        end
    end
    
    %%% mean filter
    if ~isequal(p.Results.meanFilterSize, -1)
        for angleCount = 1:angleNum
            timeSequence_thisAngle = proj_all(:,:,angleCount:angleNum:end);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/beforeMedianFilter.tiff'],optionsForTiffSave);
            end
            timeSequence_thisAngle = preprocessing_denoise_meanFilter3(timeSequence_thisAngle,p.Results.meanFilterSize);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/afterMedianFilter.tiff'],optionsForTiffSave);
            end
            proj_all(:,:,angleCount:angleNum:end) = timeSequence_thisAngle;
        end
    end
    
    %%% gaussian filter
    if ~isequal(p.Results.gaussFilterSize, -1)
        for angleCount = 1:angleNum
            timeSequence_thisAngle = proj_all(:,:,angleCount:angleNum:end);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/beforeMedianFilter.tiff'],optionsForTiffSave);
            end
            timeSequence_thisAngle = preprocessing_denoise_gaussianFilter3(timeSequence_thisAngle,p.Results.gaussFilterSize);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/afterMedianFilter.tiff'],optionsForTiffSave);
            end
            proj_all(:,:,angleCount:angleNum:end) = timeSequence_thisAngle;
        end
    end

    %%% backSubtract or not --> optionsOfBackSubtract
    if ~isequal(p.Results.backSubtract, -1)
        backSubtrackOptionCell = cell(angleNum,1);
        switch p.Results.backSubtract.method
            case {'avgSmallest','avgSmallestAbove0','mannual','someFramesAvgSmallestAbove0'}
                for angleCount = 1:angleNum
                    backSubtrackOptionCell{angleCount} = p.Results.backSubtract;
                end

            case 'set'
                for angleCount = 1:angleNum
                    backSubtrackOptionCell{angleCount} = p.Results.backSubtract;
                    backSubtrackOptionCell{angleCount}.noiseValue = p.Results.backSubtract.backLevelCell{angleCount};
                end
                
            otherwise
                error('Unknown noise determination method.');
        end
        
        for angleCount = 1:angleNum
            timeSequence_thisAngle = proj_all(:,:,angleCount:angleNum:end);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/beforeBackSubtract.tiff'],optionsForTiffSave);
            end
            timeSequence_thisAngle = preprocessing_noiseSubtraction(timeSequence_thisAngle,backSubtrackOptionCell{angleCount});
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/afterBackSubtract.tiff'],optionsForTiffSave);
            end
            proj_all(:,:,angleCount:angleNum:end) = timeSequence_thisAngle;
        end
    end

    %%% upsampling
    if ~isequal(p.Results.upsample2, -1)
        proj_all = preprocessing_stackResize(proj_all,[p.Results.upsample2,p.Results.upsample2]);
    end
    
    %%% groupMean
    if ~isequal(p.Results.groupMean, -1)
        groupSize = p.Results.groupMean;
        conv_sum_old = proj_all;
        proj_all = zeros(size(proj_all,1),size(proj_all,2),floor(size(proj_all,3)/angleNum/groupSize)*angleNum,'like',conv_sum_old);
        for angleCount = 1:angleNum
            timeSequence_thisAngle = conv_sum_old(:,:,angleCount:angleNum:end);
            timeSequence_thisAngle = groupMean(timeSequence_thisAngle,groupSize,1,3);
            if saveOn ~= 0
                saveastiff(timeSequence_thisAngle,[folderName_preprocessing{angleCount},'/afterGroupMean.tiff'],optionsForTiffSave);
            end
            proj_all(:,:,angleCount:angleNum:end) = timeSequence_thisAngle;
        end
    end
    
    %% save
    resultsSavePre = ['preprocess_',timeNow];
%     if saveOn ~= 0
        saveastiff(proj_all,[resultSave_path,'/',resultsSavePre,'.tiff']);
%         saveastiff(proj_all(:,:,1:min(angleNum,size(proj_all,3))),[resultSave_path,'/',resultsSavePre,'.tiff']);
        save([resultSave_path,'/',resultsSavePre,'.mat'],'proj_all','-v7.3');
%     end
    prepParameters = p.Results;
    save([resultSave_path,'/',resultsSavePre,'_parser.mat'],'prepParameters','returnParameters');
end
end
