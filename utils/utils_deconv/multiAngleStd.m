% Calculate time-std images for different angles
%%ELi, 20210721
function proj_all_std = multiAngleStd(proj_all,angleNum)
[w,h,a] = size(proj_all);
if angleNum == a
    proj_all_std = proj_all;
else
    proj_all_std = zeros(w,h,angleNum,'like',proj_all);
    for angleCount = 1:angleNum
        proj_all_std(:,:,angleCount) = std(proj_all(:,:,angleCount:angleNum:end),1,3);
    end
end