% rotate coordinates to a particular angle
% ELi, 20201113
function [x_coordinate_afterRotation, y_coordinate_afterRotation] = coordinateRotate(x_coordinate, y_coordinate, rotateAngle_degree)
if nargin == 2
    rotateAngle_degree = 45; %rotate 45 degrees as default
end
rotateAngle_rad = rotateAngle_degree/180*pi;
rotateMatrix = [cos(rotateAngle_rad),-sin(rotateAngle_rad);sin(rotateAngle_rad),cos(rotateAngle_rad)];

if ~isequal(size(x_coordinate),size(y_coordinate))
    error('');
end

x_coordinate_afterRotation = zeros(size(x_coordinate));
y_coordinate_afterRotation = zeros(size(y_coordinate));

for index = 1:length(x_coordinate(:))
    xy_coordinate_afterRotation = rotateMatrix*[x_coordinate(index);y_coordinate(index)];
    x_coordinate_afterRotation(index) = xy_coordinate_afterRotation(1);
    y_coordinate_afterRotation(index) = xy_coordinate_afterRotation(2);
end