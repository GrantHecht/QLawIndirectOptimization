clear; close all; clc

% Read in data
dataDir = "./data/";
kep     = readmatrix(dataDir + "kep.txt");
mee     = readmatrix(dataDir + "mee.txt");
cart    = readmatrix(dataDir + "cart.txt");
coast   = readmatrix(dataDir + "coast.txt");
time    = readmatrix(dataDir + "time.txt");

% Get last n without NaN
n = 1;
while ~isnan(time(n)) && n ~= length(time)
    n = n + 1;
end
n = n - 1;

% Strip out nans
kep     = kep(1:n,:);
mee     = mee(1:n,:);
cart    = cart(1:n,:);
coast   = coast(1:n,:);
time    = time(1:n,:);

% Get thrust and coasting arcs
cart_t = nan(size(cart));
cart_c = nan(size(cart));
sw     = false;
for i = 1:n
    if coast(i) == 0.0
        cart_t(i,:) = cart(i,:);
    else
        cart_c(i,:) = cart(i,:);
    end
    if i ~= 1
        if coast(i - 1) ~= coast(i)
            if coast(i) == 0.0
                cart_c(i,:) = cart(i,:);
            else
                cart_t(i,:) = cart(i,:);
            end
        end
    end
end

figure()
plot3(cart_t(:,1),cart_t(:,2),cart_t(:,3),"r")
hold on
plot3(cart_c(:,1),cart_c(:,2),cart_c(:,3),"b")

axis equal
grid on