% Remove the defocus item in shiftMap
% 
% ELi, 20230131
function shiftMap = removeDefocusItem(shiftMap,angDistx_normed,angDisty_normed,CAIndex)
% INPUT
%     shiftMap        - original shiftMap
%     angDistx_normed - the normalized angle distribution on pupil phase, x
%     angDisty_normed - the normalized angle distribution on pupil phase, y
%     CAIndex         - the index of centre angle

%% fit the defocus item
k1 = angDistx_normed'.*squeeze(shiftMap(2,:))+angDisty_normed'.*squeeze(shiftMap(1,:));
k2 = angDisty_normed'.*angDisty_normed'+angDistx_normed'.*angDistx_normed';
k=sum(k1(:))/sum(k2(:));
shiftMap=shiftMap-k*[angDisty_normed,angDistx_normed]';

%% set centre angle shift to 0
if exist("CAIndex",'var')
    shiftMap = shiftMap - shiftMap(:,CAIndex);
end