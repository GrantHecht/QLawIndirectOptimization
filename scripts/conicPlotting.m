clear; close all; clc

% Compute cone scale
a   = tand(90 - 30);

% Define cone parameters
[X,Y] = meshgrid(-2:0.1:2);
Z = a*sqrt(X.^2 + Y.^2);

surf(X,Y,Z)
axis equal