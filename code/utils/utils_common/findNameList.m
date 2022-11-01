% List file names in "path" that starts with "startString" and ends with "endString".
%%ELiiiiiii, 20210707
function nameList = findNameList(path,startString,endString)
%% default
if nargin == 1
    startString = [];
    endString = [];
elseif nargin == 2
    endString = [];
end
%%
oldPath = pwd;
cd(path);
result = dir([startString,'*',endString]);
nameList = {result.name};
cd(oldPath);
