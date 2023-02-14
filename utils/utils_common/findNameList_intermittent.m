% A simple function to list the files in given path
%
% NOTES
%     List file in "path" that contains intermittent strings in their names.
%     e.g. findNameList_intermittent(path,'string1','string2','string3')
%     returns file names which have the format "...string1...string2...string3..."
%
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
