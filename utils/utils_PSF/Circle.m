function result=Circle(size,position,radius)
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