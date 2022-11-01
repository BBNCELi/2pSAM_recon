function [rouOut,thetaOut]=bighopavoided_triangleModeInPolarSystem(rouNumber,rouMax,thetaNumberMin,thetaNumberMax)
%rou: linspace(0,rouMax,rouNumber+1)
%theta: linspace(0,2*pi,thetaNumber) 
% where 
% thetaNumber = thetaNumberMin(rouNumberNow==1)
% thetaNumber = thetaNumberMax(rouNumberNow==touNumber)

if (sum(size(rouNumber)~=[1 1])~=0) || (sum(size(thetaNumberMax)~=[1 1])~=0)
    error('Input error during path generation');
end

if (rouNumber<=0) || (rouNumber-fix(rouNumber)~=0) || (thetaNumberMax<=0) || (thetaNumberMax-fix(thetaNumberMax)~=0)
    error('Input error during path generation');
end

if thetaNumberMax<3
    error('');
end

rouArray = linspace(0,rouMax,rouNumber+1);
rouArray(1)=[];

%% step 0
rouOut=0;
thetaOut=0;

%% step 1
%  ---------
rouOutTemp=linspace(0,rouMax,rouNumber+1);
rouOutTemp(1)=[];
rouOut=[rouOut,rouOutTemp];
thetaOut=[thetaOut,zeros(1,rouNumber)];

%% step 2
%         |
%         |
%         |
%         |
%         |
%         |
%         |
%         |
% ---------
thetaArray = linspace(0,2*pi,thetaNumberMax+1);
thetaArray(1) = [];
thetaArray(end) = [];
rouOut=[rouOut,rouMax*ones(1,thetaNumberMax-1)];
thetaOut=[thetaOut,thetaArray];

%% step 3
%         |
%       |||
%       |||
%       |||
%       |||
%       |||
%       |||
%       |||
% ---------

for rouNumberNow = rouNumber-1:-1:1
    rouNow = rouArray(rouNumberNow);
    thetaNumberNow = getThetaNumberAccordingToRouNow(rouNow,rouMax,rouArray(1),thetaNumberMax,thetaNumberMin);
    rouOut = [rouOut,ones(1,thetaNumberNow-1)*rouNow];
    thetaArray = linspace(0,2*pi,thetaNumberNow+1);
    thetaArray(end)=[];
    thetaArray(1)=[];
    if mod(rouNumber-rouNumberNow,2)==1
        thetaOut = [thetaOut,flip(thetaArray)];
    else
        thetaOut = [thetaOut,thetaArray];
    end
end

if thetaNumberMin == 0
    rouOut(1) = [];
    thetaOut(1) = [];
end

end


function thetaNumberNow = getThetaNumberAccordingToRouNow(rouNow,rouMax,rouMin,thetaNumberMax,thetaNumberMin)
slope = (thetaNumberMax-thetaNumberMin)/(rouMax-rouMin);
thetaNumberNow = round(thetaNumberMin+slope*(rouNow-rouMin));
end