% Generate angle distribution according to given parameters
% ELi, 20220530
function [angDistx_normed, angDisty_normed] = angDistGenerator(PSFParameters,angDist_path)
angDistSaveName = angDistNameGenerator(PSFParameters.angDistParameters);

%% load existed
if exist([angDist_path,'//',angDistSaveName],'file')
    load([angDist_path,'/',angDistSaveName],'angDistx_normed','angDisty_normed');
    return;
end

%% generate new distributions
% % load reference
load([angDist_path,'//ref.mat'],...
    'reference_thetaDegree_x',...
    'reference_thetaDegree_y',...
    'reference_voltageSentToPI_X',...
    'reference_voltageSentToPI_Y');
% % rotate reference_thetaDegree_xy for some angle so that the
% % reference_thetaDegree_xy and be linear to voltageSentToPI_xy
reference_thetaDegree_y = -reference_thetaDegree_y;
tand_reference_thetaDegree_x = tand(reference_thetaDegree_x);
tand_reference_thetaDegree_y = tand(reference_thetaDegree_y);
% %-----------------try to find the best angle for rotation----------------------
% for degree = 1:360
% [reference_thetaDegreeAfterRotation_X,reference_thetaDegreeAfterRotation_Y] = coordinateRotate(reference_thetaDegree_x, reference_thetaDegree_y, degree);
% [tand_reference_thetaDegreeAfterRotation_X,tand_reference_thetaDegreeAfterRotation_Y] = coordinateRotate(tand_reference_thetaDegree_x, tand_reference_thetaDegree_y, degree);
% reference_thetaDegreeAfterRotation_X = atand(tand_reference_thetaDegreeAfterRotation_X);
% reference_thetaDegreeAfterRotation_Y = atand(tand_reference_thetaDegreeAfterRotation_Y);
% 
% [py,sy] = polyfit(reference_voltageSentToPI_Y,reference_thetaDegreeAfterRotation_Y,1);
% [px,sx] = polyfit(reference_voltageSentToPI_X,reference_thetaDegreeAfterRotation_X,1);
% s(degree) = sy.normr+sx.normr;
% end
% %-------------------------------------------
degree = 45;
[tand_reference_thetaDegreeAfterRotation_X,tand_reference_thetaDegreeAfterRotation_Y] = coordinateRotate(tand_reference_thetaDegree_x, tand_reference_thetaDegree_y, degree);
reference_thetaDegreeAfterRotation_X = atand(tand_reference_thetaDegreeAfterRotation_X);
reference_thetaDegreeAfterRotation_Y = atand(tand_reference_thetaDegreeAfterRotation_Y);
% % EXCEPT ROTATION, a stretch for voltage in one axis is also necessary.
scale = sqrt(2);
[reference_voltageSentToPI_rotated_beforeXStretch_X,reference_voltageSentToPI_rotated_beforeXStretch_Y] = coordinateRotate(reference_voltageSentToPI_X, reference_voltageSentToPI_Y, -degree);
reference_voltageSentToPI_rotated_afterXStretch_X = reference_voltageSentToPI_rotated_beforeXStretch_X*scale;
reference_voltageSentToPI_rotated_afterXStretch_Y = reference_voltageSentToPI_rotated_beforeXStretch_Y;
[reference_voltageSentToPI_afterXStretch_X,reference_voltageSentToPI_afterXStretch_Y] = coordinateRotate(reference_voltageSentToPI_rotated_afterXStretch_X, reference_voltageSentToPI_rotated_afterXStretch_Y, degree);

% % linear fit
[px,sx] = polyfit(reference_voltageSentToPI_afterXStretch_X,reference_thetaDegreeAfterRotation_X,1);
[py,sy] = polyfit(reference_voltageSentToPI_afterXStretch_Y,reference_thetaDegreeAfterRotation_Y,1);

% % switch and generate degrees
switch PSFParameters.angDistParameters.scanningMode
    case 1 %'square'
        scanningModeNow = 1;
        obj.SamplingNumber_X = PSFParameters.angDistParameters.SamplingNumber_X;
        obj.SamplingNumber_Y = PSFParameters.angDistParameters.SamplingNumber_Y;
        Interval_Y_Upper = PSFParameters.angDistParameters.Interval_Y_Upper;
        Interval_Y_Lower = -Interval_Y_Upper;
        Interval_X_Upper = PSFParameters.angDistParameters.Interval_X_Upper;
        Interval_X_Lower = -Interval_X_Upper;
        bias_X = PSFParameters.angDistParameters.bias_X;
        bias_Y = PSFParameters.angDistParameters.bias_Y;
        obj.Interval_X_Lower = Interval_X_Lower+bias_X;
        obj.Interval_X_Upper = Interval_X_Upper+bias_X;
        obj.Interval_Y_Lower = Interval_Y_Lower+bias_Y;
        obj.Interval_Y_Upper = Interval_Y_Upper+bias_Y;
    case 2 %'circ'
        scanningModeNow = 2;
        obj.rouNumber = PSFParameters.angDistParameters.lowNA_rouNumber;
        obj.thetaNumberMin = PSFParameters.angDistParameters.lowNA_thetaNumberMin;
        obj.thetaNumberMax = PSFParameters.angDistParameters.lowNA_thetaNumberMax;
        Interval_Y_Upper = PSFParameters.angDistParameters.Interval_Y_Upper;
        Interval_Y_Lower = -Interval_Y_Upper;
        Interval_X_Upper = PSFParameters.angDistParameters.Interval_X_Upper;
        Interval_X_Lower = -Interval_X_Upper;
        bias_X = PSFParameters.angDistParameters.bias_X;
        bias_Y = PSFParameters.angDistParameters.bias_Y;
        obj.Interval_X_Lower = Interval_X_Lower+bias_X;
        obj.Interval_X_Upper = Interval_X_Upper+bias_X;
        obj.Interval_Y_Lower = Interval_Y_Lower+bias_Y;
        obj.Interval_Y_Upper = Interval_Y_Upper+bias_Y;
    case 3 %'sunflower'
        scanningModeNow = 3;
        obj.NANums = PSFParameters.angDistParameters.NANums;
        Interval_Y_Upper = PSFParameters.angDistParameters.Interval_Y_Upper;
        Interval_Y_Lower = -Interval_Y_Upper;
        Interval_X_Upper = PSFParameters.angDistParameters.Interval_X_Upper;
        Interval_X_Lower = -Interval_X_Upper;
        bias_X = PSFParameters.angDistParameters.bias_X;
        bias_Y = PSFParameters.angDistParameters.bias_Y;
        obj.Interval_X_Lower = Interval_X_Lower+bias_X;
        obj.Interval_X_Upper = Interval_X_Upper+bias_X;
        obj.Interval_Y_Lower = Interval_Y_Lower+bias_Y;
        obj.Interval_Y_Upper = Interval_Y_Upper+bias_Y;
end

[target_voltageSentToPI_X,target_voltageSentToPI_Y] = creatScanPoints(obj,scanningModeNow);

target_thetaDegreeAfterRotation__X = polyval(px,target_voltageSentToPI_X);
target_thetaDegreeAfterRotation__Y = polyval(py,target_voltageSentToPI_Y);

tand_target_thetaDegreeAfterRotation_X = tand(target_thetaDegreeAfterRotation__X);
tand_target_thetaDegreeAfterRotation_Y = tand(target_thetaDegreeAfterRotation__Y);
[tand_target_thetaDegree_x,tand_target_thetaDegree_y] = coordinateRotate(tand_target_thetaDegreeAfterRotation_X, tand_target_thetaDegreeAfterRotation_Y, -degree);
target_thetaDegree_x = atand(tand_target_thetaDegree_x);
target_thetaDegree_y = atand(tand_target_thetaDegree_y);
target_thetaDegree_y = -target_thetaDegree_y;

% % convert degrees to pupil plane and normalize
angDistx_normed = sind(target_thetaDegree_y)*PSFParameters.rindexSp/PSFParameters.numAper;
angDisty_normed = sind(target_thetaDegree_x)*PSFParameters.rindexSp/PSFParameters.numAper;

%% save
save([angDist_path,'/',angDistSaveName],'angDistx_normed','angDisty_normed');
end