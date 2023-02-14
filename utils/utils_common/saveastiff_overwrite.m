% A simple expansion for saveastiff that set overwriting as default and 
% save .mat at the same time
%
%%ELiiiiiii, 20211215
function saveastiff_overwrite(data, path, matSaveOn,tifSaveOn)
if ~exist('matSaveOn','var') && ~exist('tifSaveOn','var')
    options.overwrite = true;
    saveastiff(data,path,options);
elseif exist('matSaveOn','var') && ~exist('tifSaveOn','var')
    if matSaveOn ~= 0
        save([path,'.mat'],'data','-v7.3');
    end
elseif exist('matSaveOn','var') && exist('tifSaveOn','var')
    if matSaveOn ~= 0
        save([path,'.mat'],'data','-v7.3');
    end
    if tifSaveOn ~= 0
        options.overwrite = true;
        saveastiff(data,path,options);
    end
end