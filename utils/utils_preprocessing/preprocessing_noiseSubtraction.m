% An exact-value-subtraction is necessary for both noise reduction and background removal.
% The main problem is to decide the exact value under which should be considered as noise or backgound.
% As the BACKGROUND value is usually inhomogeneous among the image whcih is quite hard to remove,
% we here termed this function as NOISE subtraction,
% although a simple subtraction can't really remove the noise.

% ELiiiiiii, 20210123. 'avgSmallest' ''avgSmallestAbove0' 'mannual' 'set'
% ELiiiiiii, 20210531. Add 'someFramesAvgSmallestAbove0'
function EI_after = preprocessing_noiseSubtraction(EI_before,options)
% Inputs:
%     EI_before: a 3-D stack in which the frames are all captured under the SAME experimental conditions.
%     options     
%         method: the method to compute noise value.
%             'avgSmallest': the average of the SEVERAL smallest values. SEVERAL: a percentage given by options.percentage
%             'avgSmallestAbove0': the average of the SEVERAL smallest positive values. SEVERAL: a percentage given by options.percentage
%             'mannual': input the noise position according to plotted intensity distribution.
%             'set': using given noise value given by options.noiseValue.
% Outputs:
%     EI_after: a 3-D stack with the same size of EI_before.

%% dimension judgement
dim_EI = ndims(EI_before);
switch(dim_EI)
    case 1
        warning('SUSPICIOUS: Subtracting noise on an array...');
    case 2
        disp('Subtracting noise on the image!');
    case 3
        disp('Subtracting noise on the whole stack!');
    otherwise
        error('Noise subtraction for 4-D or higher-dimension input is temporarily unsupported.');
end

%% preparations
totalPixelsNum = length(EI_before(:));

%% calcute noise level
switch(options.method)
    case 'avgSmallest'
        portion = options.percentage;%Set the average of SEVERAL smallest values as noise.
        noisePixelsNum = round(portion * totalPixelsNum);
        noise = mean( mink(EI_before(:),noisePixelsNum,'ComparisonMethod','auto') );
        
    case 'avgSmallestAbove0'
        portion = options.percentage;%Set the average of SEVERAL smallest positive values as noise.
        EI_before_array = EI_before(:);
        EI_before_array(EI_before_array<0) = [];
        noisePixelsNum = round(portion * totalPixelsNum);
%         noise = mean( mink(EI_before_array(:),noisePixelsNum,'ComparisonMethod','auto') );
        noise = mean( mink(EI_before_array(:),noisePixelsNum) );
        
    case 'mannual'
%         EI_before_array = mink(EI_before(:),totalPixelsNum,'ComparisonMethod','auto');
        EI_before_array = mink(EI_before(:),totalPixelsNum);
        figure;plot(EI_before_array);
        noisePixelsNum = input('Please input the position of the first inflection: ');
        noise = EI_before_array(noisePixelsNum);
        
    case 'set'
        noise = options.noiseValue;
        
    case 'someFramesAvgSmallestAbove0'
        portion = options.percentage;%Set the average of SEVERAL smallest positive values as noise.
        someFramesToCountAvg = options.someFramesIndex;
        EI_before_array = EI_before(:,:,someFramesToCountAvg);
        EI_before_array = EI_before_array(:);
        EI_before_array(EI_before_array<0) = [];
        noisePixelsNum = round(portion * totalPixelsNum);
%         noise = mean( mink(EI_before_array(:),noisePixelsNum,'ComparisonMethod','auto') );
        noise = mean( mink(EI_before_array(:),noisePixelsNum) );
        
    otherwise
        error('Unknown noise determination method.');
        
end

%% subtraction
disp(['Noise level: ', num2str(noise)]);
EI_after = EI_before - noise;
EI_after(EI_after<0) = 0;

close all;
