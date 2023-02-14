% This function registers the input stack slice-by-slice to fixedFrame
% 
% ELiiiiiii, 20210317.
% ELiiiiiii, 20210324, add inputTemplate
function [outputStack,template] = preprocessing_2DRegistration(inputStack,inputTemplate)
% Inputs:
%     inputStack: 3-D stack to be registered.
%     inputTemplate: the template that is used for motion correction.
% Outputs:
%     outputStack: 3-D stack after registration.
%
%%% add normcorre functions first
%%% https://github.com/flatironinstitute/NoRMCorre

%% perform motion correction using normcorre
options_rigid = NoRMCorreSetParms('d1',size(inputStack,1),'d2',size(inputStack,2),'bin_width',100,'max_shift',50,'us_fac',50,'init_batch',20);
if nargin < 2 || sum(inputTemplate(:))==0 %Obtain template automatically
    [outputStack,shifts,template,options,col_shift] = normcorre(inputStack,options_rigid);
else
    [outputStack,shifts,template,options,col_shift] = normcorre(inputStack,options_rigid,inputTemplate);
end
