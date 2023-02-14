% Calculate time-std images at different angles
% 
% ELi, 20210721
% ELi, 20230207, add comments
function projs_std = multiAngleStd(projs,angleNum)
% INPUT
%     projs     - time-series projections
%     angleNum  - number of angles contained in proj_all

%% preparations
[w,h,a] = size(projs);

%% time-std angle by angle
if angleNum == a % if no need to calculate std, return proj_all directly
    projs_std = projs;
else
    projs_std = zeros(w,h,angleNum,'like',projs);
    for angleCount = 1:angleNum
        projs_std(:,:,angleCount) = std(projs(:,:,angleCount:angleNum:end),1,3);
    end
end