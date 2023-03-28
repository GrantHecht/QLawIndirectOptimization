clear; close all; clc

addpath("./scripts/utils/");

% Data name
%dn      = "cp";
dn      = "cdp";
%dn      = "ctd";

% Figure folder
figFolder = "C:\Users\grant\Documents\Projects\Astro-2023-Low-Thrust-Q-Law-Control-With-Sun-Angle-Constraint\Figures\";

% Read in data
dataDir = "./data/";
kep     = readmatrix(dataDir + "kep_" + dn + ".txt");
mee     = readmatrix(dataDir + "mee_" + dn + ".txt");
cart    = readmatrix(dataDir + "cart_" + dn + ".txt");
coast   = readmatrix(dataDir + "coast_" + dn + ".txt");
eclipse = readmatrix(dataDir + "eclipse_" + dn + ".txt");
angles  = readmatrix(dataDir + "angles_" + dn + ".txt");
thrust  = readmatrix(dataDir + "thrust_" + dn + ".txt");
time    = readmatrix(dataDir + "time_" + dn + ".txt");
kept    = readmatrix(dataDir + "kept_" + dn + ".txt");
consts  = readmatrix(dataDir + "consts_" + dn + ".txt");
sunangs = readmatrix(dataDir + "sunangles_" + dn + ".txt");

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

% Print info
fprintf("TOF:    %1.4f days\n", time(end));
fprintf("DeltaM: %1.4f kg\n", kep(1,7) - kep(end,7));


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
plotKeplerianOrbit(kept, 360, consts(1))
xlabel("x")
ylabel("y")

axis equal
grid on

figure()
tiledlayout(3,1)

% X-Y Plot
nexttile([2,1])
plot(cart_t(:,1),cart_t(:,2),"r")
hold on
plot(cart_e(:,1),cart_e(:,2),"g")
plot(cart_c(:,1),cart_c(:,2),"b")
plotKeplerianOrbit(kept, 360, consts(1))

% Setup axes
xlabel("x, km", "Interpreter", "Latex")
ylabel("y, km", "Interpreter", "Latex")
axis equal
grid on

% Setup legend
leg1 = legend("Thrusting","Ecclipsed","Coasting","Location","northwest");

% Publication plotting settings
set(gca, "fontname", "Times New Roman", "fontsize", 10)

% Y-Z Plot
nexttile
plot(cart_t(:,2),cart_t(:,3),"r")
hold on
plot(cart_e(:,2),cart_e(:,3),"g")
plot(cart_c(:,2),cart_c(:,3),"b")

% Setup axes
xlabel("y, km", "Interpreter", "Latex")
ylabel("z, km", "Interpreter", "Latex")

axis equal
grid on

% Publication plotting settings
set(gca, "fontname", "Times New Roman", "fontsize", 10)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,5,4])
set(gcf, "PaperPositionMode","Manual")

% Print figure
print(figFolder + "traj_full_" + dn + ".eps", "-depsc", "-r300")

figure()
tiledlayout(4,1);
for j = 1:4
    if j == 4
        i = 7;
    %elseif j == 5
    %    i = 7;
    else
        i = j;
    end
    nexttile
    plot(time, kep(:,i))
    if j == 1
        ylabel("a, km", "Interpreter", "Latex")
        ylim([0.0, 6.7e4])
    elseif j == 2
        ylabel("e", "Interpreter", "Latex")
        ylim([-0.1, 0.85])
        yticks([0,0.3,0.6])
    elseif j == 3
        ylabel("i, deg", "Interpreter", "Latex")
        ylim([-0.1, 0.55])
    elseif j == 4 
        xlabel("time, days", "Interpreter", "Latex")
        ylabel("m, kg", "Interpreter", "Latex")
        ylim([850, 1250])
    end
    grid on
    set(gca, "fontname", "Times New Roman", "fontsize", 10)
end

% Publication plotting settings
set(gca, "fontname", "Times New Roman", "fontsize", 10)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,5.0,3.5])
set(gcf, "PaperPositionMode","Manual")

% Print figure
print(figFolder + "timeHists_" + dn + ".eps", "-depsc", "-r300")

figure()
scatter(time,sunangs,".")
xlabel("time, days", "Interpreter", "Latex")
ylabel("sun angle, deg", "Interpreter", "Latex")

%ylim([0,180])
grid on

% Publication plotting settings
set(gca, "fontname", "Times New Roman", "fontsize", 10)
%if dn == "cp"
    set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,6.0,2.5])
%else
%    set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,3.0,2.5])
%end
set(gcf, "PaperPositionMode","Manual")

% Print figure
print(figFolder + "sunAngles_" + dn + ".eps", "-depsc", "-r300")
