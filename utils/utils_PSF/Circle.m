% A simple functon to generate a circle
%
% ELi, 20230213, add comments
function result=Circle(size,position,radius)
% INPUTS:
%     size    - size(result)
%     postion - centre of the circle
%     radius  - the inner and outer radius in pixel of the circle

result=zeros(size(1),size(2));
for i=1:size(1)
    for j=1:size(2)
        if (((i-position(1))^2+(j-position(2))^2)<=radius(2)^2)&&...
            (((i-position(1))^2+(j-position(2))^2)>=radius(1)^2)
        result(i,j)=1;
        end
    end
end
end