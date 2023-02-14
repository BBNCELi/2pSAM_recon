% A simple expansion for "mapminmax" that accepts stacks as input
%
% ELi, 20230214, add comments
function originalStack = mapminmax_wholeStack(originalStack, minValue, maxValue)
%% max and min
originalMax = max(originalStack(:));
originalMin = min(originalStack(:));

%% linear stretch
originalStack = (originalStack - originalMin) * ((maxValue - minValue)/(originalMax - originalMin)) + minValue;
end