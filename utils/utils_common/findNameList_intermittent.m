% List file names in "path" that contains intermittent strings in their names.
%     e.g. findNameList_intermittent(path,'string1','string2','string3')
%     will return file names which have the format
%     "...string1...string2...string3..."
%%ELiiiiiii, 20210708
function nameList = findNameList_intermittent(path,varargin)
%% set name format
nameFormat = '*';
for i = 1:length(varargin)
    nameFormat = [nameFormat,varargin{i},'*'];
end

%%
oldPath = pwd;
cd(path);
result = dir(nameFormat);
nameList = {result.name};
cd(oldPath);
