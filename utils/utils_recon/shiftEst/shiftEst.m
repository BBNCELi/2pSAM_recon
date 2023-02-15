% Calculate the shift map and then the disparity map between imgRef and imgMov,
% according to given patch options.
% 
% NOTES
%     When calculating the disparity map from the shift map, this function:
%       uses linear interpolation for the overlapping areas
%       uses nearest extrapolation for the sidelobe
% 
% ELi, 20230207
function [shiftMap, dispMap, imgMov_corrected] = shiftEst(imgRef, imgMov, opts)
% INPUT
%     imgRef    - reference image
%     imgMov    - translated image
%     opts      - options for shift estimation
%                 .shiftEstFunc: function used for shift estimation
%                 .patchN: number of patches for multi-site shift estimation
%                 .patchOvFactor: overlap factor for multi-site shift estimation
%                 .minPatchSize: minimum patch size for multi-site shift estimation
%                 .sidelobe: sidelobe in pixel for shift estimation

%% check inputs; extract opts for better use; other preparations
if ~exist('opts','var'); opts = struct; end
if ~isfield(opts,'shiftEstFunc'); opts.shiftEstFunc = 'imregtform'; end
if ~isfield(opts,'patchN'); opts.patchN = 1; end
if ~isfield(opts,'patchOvFactor'); opts.patchOvFactor = 0; end
if ~isfield(opts,'minPatchSize'); opts.minPatchSize = 1; end
if ~isfield(opts,'sidelobe'); opts.sidelobe = 0; end

shiftEstFunc = opts.shiftEstFunc;
patchN = opts.patchN;
patchOvFactor = opts.patchOvFactor;
minPatchSize = opts.minPatchSize;
sidelobe = opts.sidelobe;

if ~ismatrix(imgRef)
    error('this function supports only images as inputs');
end
if ~isequal(size(imgRef), size(imgMov))
    error('imgRef and imgMov should have the same size');
end
[img_r,img_c] = size(imgRef);

%% patch construction and shift estimation
[imgRef_patch,patchParas] = mat2patch(imgRef,patchN,patchOvFactor,minPatchSize,sidelobe);
[imgMov_patch,patchParas] = mat2patch(imgMov,patchN,patchOvFactor,minPatchSize,sidelobe);
patchN_x = length(patchParas.xx_s_ov_sidelobe);
patchN_y = length(patchParas.yy_s_ov_sidelobe);
shiftMap = zeros(patchN_x,patchN_y,2); % x&&y

for i = 1 : patchN_x
    for j = 1 : patchN_y
        img1 = imgRef_patch{i,j};
        img2 = imgMov_patch{i,j};

        %%%estimate shifts using shiftEstFunc to avoid using switch
%         [shiftMap(i,j,1),shiftMap(i,j,2)] = ...
%             eval(['shiftEst_',shiftEstFunc,'(img1,img2)']);
        switch shiftEstFunc
            case 'corr'
                [shiftMap(i,j,1),shiftMap(i,j,2)] = ...
                shiftEst_corr(img1,img2);
            case 'Foroosh2002'
                [shiftMap(i,j,1),shiftMap(i,j,2)] = ...
                shiftEst_Foroosh2002(img1,img2);
            case 'imregtform'
                [shiftMap(i,j,1),shiftMap(i,j,2)] = ...
                shiftEst_imregtform(img1,img2);
            case 'phaseCorr'
                [shiftMap(i,j,1),shiftMap(i,j,2)] = ...
                shiftEst_phaseCorr(img1,img2);
            otherwise
                error('Unrecognized shift estimation function');
        end
    end
end

%% interpolate the shift map to get disparity map if asked
if nargout > 1
    %%%initialize
    dispMap = nan(img_r,img_c,2); % x&&y
    if patchN_x == 1 && patchN_y == 1
        %%%%when patchN == 1, no need to interpolate
        dispMap(:,:,1) = shiftMap(:,:,1);
        dispMap(:,:,2) = shiftMap(:,:,2);
    else
        %%%patch paras
        xx_s = patchParas.xx_s_ov_sidelobe;
        xx_f = patchParas.xx_f_ov_sidelobe;
        yy_s = patchParas.yy_s_ov_sidelobe;
        yy_f = patchParas.yy_f_ov_sidelobe;
    
        %%%weight_patchx
        weight_patchx = zeros(img_r,patchN_x);
        %%%%for no-overlapping areas, weight_patchx = 1
        for i = 1:img_r
            isInPatch_x = (i>=xx_s) .* (i<=xx_f);
            numInPatch_x = sum(isInPatch_x);
            if numInPatch_x == 1
                patch_x = isInPatch_x == 1;
                weight_patchx(i,patch_x) = 1;
            end
        end
        %%%%for overlapping areas, linear interpolation
        for patchCount = 2:length(xx_s)
            sect = xx_s(patchCount) : xx_f(patchCount-1);
            sect_l = length(sect);
            weight_patchx(sect,patchCount-1) = linspace(1,0,sect_l);
            weight_patchx(sect,patchCount) = linspace(0,1,sect_l);
        end
        %%%weight_patchy
        %%%%for no-overlapping areas, weight_patchy = 1
        weight_patchy = zeros(img_c,patchN_y);
        for j = 1:img_c
            isInPatch_y = (j>=yy_s) .* (j<=yy_f);
            numInPatch_y = sum(isInPatch_y);
            if numInPatch_y == 1
                patch_y = isInPatch_y == 1;
                weight_patchy(j,patch_y) = 1;
            end
        end
        %%%%for overlapping areas, linear interpolation
        for patchCount = 2:length(yy_s)
            sect = yy_s(patchCount) : yy_f(patchCount-1);
            sect_l = length(sect);
            weight_patchy(sect,patchCount-1) = linspace(1,0,sect_l);
            weight_patchy(sect,patchCount) = linspace(0,1,sect_l);
        end
        %%%start calculating disparity map
        for i = xx_s(1) : xx_f(end) %skip sidelobe at this moment
            for j = yy_s(1) : yy_f(end) %skip sidelobe at this moment
                weight_x = weight_patchx(i,:);
                weight_y = weight_patchy(j,:);
                weight = weight_x' * weight_y;
                dispMap(i,j,1) = sum(shiftMap(:,:,1) .* weight,'all');
                dispMap(i,j,2) = sum(shiftMap(:,:,2) .* weight,'all');
            end
        end
    
        %%%for sidelobes, use nearest interpolation
        dispMap = fillmissing(dispMap,'nearest');
        dispMap(:,:,1) = fillmissing(dispMap(:,:,1)','nearest')';
        dispMap(:,:,2) = fillmissing(dispMap(:,:,2)','nearest')';
    end
end

%% shift imgMov according to the estimated dispMap if asked
if nargout > 2
    imgMov_corrected = imtranslate_disparityMap(imgMov, dispMap);
end