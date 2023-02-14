% Estimate shift between projections and centre angle projection
% 
% ELi, 20230131
function shiftMap = shiftToCA(projs,CAIndex,shiftEstFunc)
% INPUT
%     projs           - projections from all angles
%     CAIndex         - the index of centre angle
%     shiftEstFunc    - shift estimator

%% check input and size
[proj_r,proj_c,angleNum] = size(projs);
shiftMap = zeros(2, angleNum);
if ~exist("shiftEstFunc",'var'); shiftEstFunc = 'corr'; end

%% estimate shifts one by one
proj_CA = projs(:,:,CAIndex);
for angleNow = 1:angleNum
    projNow = projs(:,:,angleNow);
    [shiftMap(1,angleNow),shiftMap(2,angleNow)] = ...
        eval(['shiftEst_',shiftEstFunc,'(projNow,proj_CA)']);
end