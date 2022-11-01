function originalStack = mapminmax_wholeStack(originalStack, minValue, maxValue)
originalMax = max(originalStack(:));
originalMin = min(originalStack(:));

originalStack = (originalStack - originalMin) * ((maxValue - minValue)/(originalMax - originalMin)) + minValue;
end