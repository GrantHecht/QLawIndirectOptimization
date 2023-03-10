clear; close all; clc

addpath("./scripts/utils/");

% Read in data
dataDir = "./data/";
kep     = readmatrix(dataDir + "kep.txt");
mee     = readmatrix(dataDir + "mee.txt");
cart    = readmatrix(dataDir + "cart.txt");
coast   = readmatrix(dataDir + "coast.txt");
eclipse = readmatrix(dataDir + "eclipse.txt");
angles  = readmatrix(dataDir + "angles.txt");
thrust  = readmatrix(dataDir + "thrust.txt");
time    = readmatrix(dataDir + "time.txt");
kept    = readmatrix(dataDir + "kept.txt");
consts  = readmatrix(dataDir + "consts.txt");
sunangs = readmatrix(dataDir + "sunangles.txt");

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
sunangs = sunangs(1:n,:);

% Get thrust and coasting arcs
cart_t = nan(size(cart));
cart_e = nan(size(cart));
cart_c = nan(size(cart));
kep_t  = nan(size(kep));
kep_e  = nan(size(kep));
kep_c  = nan(size(kep));
sw     = false;
for i = 1:n
    tnow = false;
    enow = false;
    cnow = false;
    if coast(i) == 0.0 && eclipse(i) == 0.0
        cart_t(i,:) = cart(i,:);
        kep_t(i,:) = kep(i,:);
        tnow = true;
    elseif eclipse(i) ~= 0.0 
        cart_e(i,:) = cart(i,:);
        kep_e(i,:) = kep(i,:);
        enow = true;
    else
        cart_c(i,:) = cart(i,:);
        kep_c(i,:) = kep(i,:);
        cnow = true;
    end
    if i ~= 1
        tlast = false;
        elast = false;
        clast = false;
        if coast(i - 1) == 0.0 && eclipse(i - 1) == 0.0
            tlast = true;
        elseif eclipse(i - 1) ~= 0.0
            elast = true;
        else 
            clast = true;
        end
        if tnow == false && tlast == true
            cart_t(i,:) = cart(i,:);
            kep_t(i,:) = kep(i,:);
        elseif enow == false && elast == true
            cart_e(i,:) = cart(i,:);
            kep_e(i,:) = kep(i,:);
        elseif cnow == false && clast == true 
            cart_c(i,:) = cart(i,:);
            kep_c(i,:) = kep(i,:);
        end
    end
end

figure()
plot3(cart_t(:,1),cart_t(:,2),cart_t(:,3),"r")
hold on
plot3(cart_c(:,1),cart_c(:,2),cart_c(:,3),"b")
plot3(cart_e(:,1),cart_e(:,2),cart_e(:,3),"g")
scatter3(cart([1,end],1),cart([1,end],2),cart([1,end],3));
plotKeplerianOrbit(kept, 360, consts(1))
xlabel("x")
ylabel("y")

axis equal
grid on

figure()
tiledlayout(5,1);
for j = 1:5
    if j == 4
        i = 6;
    elseif j == 5
        i = 7;
    else
        i = j;
    end
    nexttile
    scatter(time, kep_t(:,i),1.5,"r","filled")
    hold on
    scatter(time, kep_c(:,i),1.5,"b","filled")
    scatter(time, kep_e(:,i),1.5,"g","filled")
    hold off
end

figure()
scatter(time,sunangs,".")
