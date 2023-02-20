clear; close all; clc
a   = tand(90 - 20);

[X,Y] = meshgrid(-2:0.1:2);
Z = sqrt(X.^2 + Y.^2/4)/a;
surf(X,Y,Z)
view(8,2)
xlabel 'x'
ylabel 'y'
zlabel 'z'
axis equal
