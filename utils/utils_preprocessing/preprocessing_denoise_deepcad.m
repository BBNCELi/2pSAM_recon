%Use deepcad to denoise input image sequence
%%ELiiiiiii, 20210625
% 
%Modification for better use
%%ELiiiiiii, 20211029
%%ELiiiiiii, 20220302, update Deepcad-RT
function outputStack = preprocessing_denoise_deepcad(inputStack,model)
% INPUTS:
%     inputStack: 3D stack
%     model:deepcad model that will be used.
%            For example, if model = 'cwin', then this function will use model "cwin.engine" to denoise inputStack.
%            NOTICE: The system environment should be correctly set to make sure that "cwin.cpp" runs fine.
% 
%%% add deepcad-rt functions first
%%% https://github.com/cabooster/DeepCAD-RT

%%
% inputStack(inputStack<0) = 0;
mean_inputStack = mean(mean(inputStack,2),1);
%%%
inputStack = double(inputStack);%Because .mexw64 file do not accept double input
[hei,wid,frame_input] = size(inputStack);
INPUT_SEQ_H = 200; % Patch size in height
INPUT_SEQ_W = 200; % Patch size in width
INPUT_SEQ_S = 80; % Patch size in slice
ov_factor=0.4; %the overlap factor between two adjacent patches
GAP_H = round(INPUT_SEQ_H*(1-ov_factor));% Patch gap in height
GAP_W = round(INPUT_SEQ_W*(1-ov_factor));% Patch gap in width
GAP_S = round(INPUT_SEQ_S*(1-ov_factor));% Patch gap in slice
normalize_factor=1;
model_fname = which([model,'.engine']);
if isempty(model_fname)
    error(['No such model for deepcad: ',model,'.engine']);
end
outputStack = deepcad_core_saveinmatlab(inputStack,wid,hei,frame_input,INPUT_SEQ_H, INPUT_SEQ_W, INPUT_SEQ_S, GAP_H, GAP_W, GAP_S,model_fname,normalize_factor);
outputStack = reshape(outputStack,hei,wid,frame_input);
%%%
% outputStack(outputStack<0) = 0;
mean_outputStack = mean(mean(outputStack,2),1);
% outputStack = mapminmax_wholeStack(outputStack,min(inputStack(:)),max(inputStack(:)));
outputStack = outputStack./mean_outputStack.*mean_inputStack;
outputStack = single(outputStack);%Change back to single
