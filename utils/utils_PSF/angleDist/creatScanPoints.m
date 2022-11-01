% This function was part of the system control unit
% Modified for PSF simulation, ELi, 20220530
function [voltageSentToPI_X,voltageSentToPI_Y] = creatScanPoints(obj,scanningModeNow)
    if scanningModeNow == 1
        Interval_X=[obj.Interval_X_Lower obj.Interval_X_Upper];         %maximum:[-10 10](both MEMS&6363)
        Interval_Y=[obj.Interval_Y_Lower obj.Interval_Y_Upper];
        SamplingPointsInTotal=obj.SamplingNumber_X*obj.SamplingNumber_Y;
        SamplingValues_X=linspace(Interval_X(1),Interval_X(2),obj.SamplingNumber_X);     %evenly distrubuted
        SamplingValues_Y=linspace(Interval_Y(1),Interval_Y(2),obj.SamplingNumber_Y);
        [SamplingPath_X,SamplingPath_Y]=bighopavoided(obj.SamplingNumber_X,obj.SamplingNumber_Y);
        SamplingPoints_X=[];
        SamplingPoints_Y=[];
        for i=1:SamplingPointsInTotal
            SamplingPoints_X(i)=SamplingValues_X(SamplingPath_X(i));
            SamplingPoints_Y(i)=SamplingValues_Y(SamplingPath_Y(i));
        end
        [ScanPosition_X_temp,ScanPosition_Y_temp]=coordinateRotate(SamplingPoints_X,SamplingPoints_Y,45);  %Rotate
    elseif scanningModeNow == 2
        lowNA_rouMax = 1;%put shift_x,shift_y in [-1,1]
        [shift_rou,shift_theta]=bighopavoided_triangleModeInPolarSystem(obj.rouNumber,lowNA_rouMax,obj.thetaNumberMin,obj.thetaNumberMax);
        for i = 1 : length(shift_rou)
                ScanPosition_X_temp(i) = shift_rou(i).*cos(shift_theta(i));
                ScanPosition_Y_temp(i) = shift_rou(i).*sin(shift_theta(i));
        end
        %linear stretching
        ScanPosition_X_temp = (ScanPosition_X_temp-(-1)).*(obj.Interval_X_Upper - obj.Interval_X_Lower)./2 + obj.Interval_X_Lower;
        ScanPosition_Y_temp = (ScanPosition_Y_temp-(-1)).*(obj.Interval_Y_Upper - obj.Interval_Y_Lower)./2 + obj.Interval_Y_Lower;
        [ScanPosition_X_temp,ScanPosition_Y_temp]=coordinateRotate(ScanPosition_X_temp,ScanPosition_Y_temp,45);  %Rotate
    elseif scanningModeNow == 3
        [ScanPosition_X_temp,ScanPosition_Y_temp] = sunflowerDistribution(obj.NANums, 0);
        %linear stretching
        ScanPosition_X_temp = (ScanPosition_X_temp-(min(ScanPosition_X_temp))).*(obj.Interval_X_Upper - obj.Interval_X_Lower)./2 + obj.Interval_X_Lower;
        ScanPosition_Y_temp = (ScanPosition_Y_temp-(min(ScanPosition_Y_temp))).*(obj.Interval_Y_Upper - obj.Interval_Y_Lower)./2 + obj.Interval_Y_Lower;
        [ScanPosition_X_temp,ScanPosition_Y_temp]=coordinateRotate(ScanPosition_X_temp,ScanPosition_Y_temp,45);  %Rotate
    end

    ScanPosition_X_withoutMultipleFrames = ScanPosition_X_temp;
    ScanPosition_Y_withoutMultipleFrames = ScanPosition_Y_temp;
    
    voltageSentToPI_X = ScanPosition_X_withoutMultipleFrames./1.75;
    voltageSentToPI_Y = ScanPosition_Y_withoutMultipleFrames./1.75;
end