clear; close all; clc

addpath("./scripts/utils/");

% Read in data
dataDir = "./data/";
kep     = readmatrix(dataDir + "kep.txt");
mee     = readmatrix(dataDir + "mee.txt");
cart    = readmatrix(dataDir + "cart.txt");
coast   = readmatrix(dataDir + "coast.txt");
angles  = readmatrix(dataDir + "angles.txt");
thrust  = readmatrix(dataDir + "thrust.txt");
time    = readmatrix(dataDir + "time.txt");
kept    = readmatrix(dataDir + "kept.txt");
consts  = readmatrix(dataDir + "consts.txt");

% Sun direction
toSun   = [1.0; 0.0; 0.0];

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
angles  = angles(1:n,:);
thrust  = thrust(1:n);
time    = time(1:n,:);

% Get thrust and coasting arcs
cart_t = nan(size(cart));
cart_c = nan(size(cart));
stimes = [];
sw     = false;
for i = 1:n
    if coast(i) == 0.0
        cart_t(i,:) = cart(i,:);
    else
        cart_c(i,:) = cart(i,:);
    end
    if i ~= 1
        if coast(i - 1) ~= coast(i)
            stimes = [stimes;time(i)];
            if coast(i) == 0.0
                cart_c(i,:) = cart(i,:);
            else
                cart_t(i,:) = cart(i,:);
            end
        end
    end
end

% Compute angle between thrust and sun
sunAngles = zeros(length(time),1);
for i = 1:length(time)
    if coast(i) == 0
        % Construct rotation from LVLH to Inertial
        ur    = cart(i,1:3) / norm(cart(i,1:3));
        uh    = cross(cart(i,1:3), cart(i,4:6));
        uh    = uh / norm(uh);
        ut    = cross(uh,ur);
        A     = [ur', ut', uh'];
    
        % Construct thrust unit vector
        ut    = A*[cos(angles(i,2))*sin(angles(i,1)); 
                   cos(angles(i,2))*cos(angles(i,1));
                   sin(angles(i,2))];
    
        dprod = toSun'*ut;
        sunAngles(i) = acosd(abs(dprod) / (norm(toSun)*norm(ut)));
    else
        sunAngles(i) = NaN;
    end
end

% figure()
% plot(mee(:,1))
% 
% figure()
% plot(mee(:,7))

figure()
plot3(cart_t(:,1),cart_t(:,2),cart_t(:,3),"r")
hold on
plot3(cart_c(:,1),cart_c(:,2),cart_c(:,3),"b")
scatter3(cart([1,end],1),cart([1,end],2),cart([1,end],3));
plotKeplerianOrbit(kept, 360, consts(1))
xlabel("x")
ylabel("y")

axis equal
grid on

figure()
scatter(time,sunAngles,".")
