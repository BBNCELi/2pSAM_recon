% Load the first variable in .mat file and return it back
% For using "load" more conveniently

%ELiiiiiii, 20211211
function var = loadFirstVariable(matFile)
    matStruct=load(matFile);
    matCell=fieldnames(matStruct);
    variableName=matCell{1};%load the first variable
    eval(['var=matStruct.',variableName,';']);%name it 'var'
end
