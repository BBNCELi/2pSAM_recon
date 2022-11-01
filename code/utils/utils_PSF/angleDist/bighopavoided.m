function [xout,yout]=bighopavoided(xin,yin)

if (sum(size(xin)~=[1 1])~=0) || (sum(size(yin)~=[1 1])~=0)
    error('Input error during path generation');
end

if (xin<=0) || (xin-fix(xin)~=0) || (yin<=0) || (yin-fix(yin)~=0)
    error('Input error during path generation');
end

%% step 1
%  ---------
xout=1:1:xin;
yout=ones(1,xin);

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
xout=[xout,xin*ones(1,yin-1)];
yout=[yout,2:yin];

%% step 3
% --------|
% --------|
% --------|
% --------|
%         |
%         |
%         |
%         |
% ---------
if mod(yin,2)==0
    for m=0:1:(yin-4)/2
        xout=[xout,xin-1:-1:1,1:1:xin-1];
        yout=[yout,(yin-2*m)*ones(1,xin-1),(yin-2*m-1)*ones(1,xin-1)];
    end
    xout=[xout,xin-1:-1:1];
    yout=[yout,2*ones(1,xin-1)];
    
elseif mod(yin,2)==1
    for m=0:1:(yin-5)/2
        xout=[xout,xin-1:-1:1,1:1:xin-1];
        yout=[yout,(yin-2*m)*ones(1,xin-1),(yin-2*m-1)*ones(1,xin-1)];
    end
   
%%step 4
% --------|
% --------|
% --------|
% --------|
% --------|
% --------|
%     |||||
%     |||||
% ---------
    if mod(xin,2)==0
        for n=0:1:(xin-4)/2
            xout=[xout,(xin-2*n-1),(xin-2*n-1),(xin-2*n-2),(xin-2*n-2)];
            yout=[yout,3,2,2,3];
        end
        xout=[xout,1,1];
        yout=[yout,3,2];
    elseif mod(xin,2)==1
        for n=0:1:(xin-5)/2
            xout=[xout,(xin-2*n-1),(xin-2*n-1),(xin-2*n-2),(xin-2*n-2)];
            yout=[yout,3,2,2,3];
        end
        
%%step 5
% --------|
% --------|
% --------|
% --------|
% --------|
% --------|
%     |||||
%     |||||
% ---------
        xout=[xout,2,1,1,2];
        yout=[yout,3,3,2,2];
    end
end