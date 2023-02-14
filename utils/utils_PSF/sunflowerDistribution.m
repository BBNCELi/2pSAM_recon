% A simple function to generate a sunflower-like distribution
%
%Eli, 20230213, add comments
function [xout,yout]=sunflowerDistribution(n, alpha)   %  example: n=500, alpha=2
    xout=[];
    yout=[];
%     hold on
    b = round(alpha*sqrt(n));      % number of boundary points
    phi = (sqrt(5)+1)/2;           % golden ratio
    for k=1:n
        r = radius(k,n,b);
        theta = 2*pi*k/phi^2;
%         plot(r*cos(theta), r*sin(theta), 'r*');
        xout = [xout,r*cos(theta)];
        yout = [yout,r*sin(theta)];
    end
end

function r = radius(k,n,b)
    if k>n-b
        r = 1;            % put on the boundary
    else
        r = sqrt(k-1/2)/sqrt(n-(b+1)/2);     % apply square root
    end
end