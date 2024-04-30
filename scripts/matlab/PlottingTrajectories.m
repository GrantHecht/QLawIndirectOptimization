clear; close all; clc

addpath("./scripts/utils/");

% Data name
%dn     = "";
dn      = "_cp";
%dn      = "_cdp";
%dn      = "_ctd";

% Figure folder
%figFolder = "C:\Users\grant\Documents\Projects\Astro-2023-Low-Thrust-Q-Law-Control-With-Sun-Angle-Constraint\Figures\";
%figFolder = "/Users/granthec/Documents/Projects/Astro 2023/Astro-2023-Low-Thrust-Q-Law-Control-With-Sun-Angle-Constraint\Figures\";
figFolder = "./plots/";

% Read in data
dataDir = "./data/";
kep     = readmatrix(dataDir + "kep" + dn + ".txt");
mee     = readmatrix(dataDir + "mee" + dn + ".txt");
cart    = readmatrix(dataDir + "cart" + dn + ".txt");
coast   = readmatrix(dataDir + "coast" + dn + ".txt");
eclipse = readmatrix(dataDir + "eclipse" + dn + ".txt");
angles  = readmatrix(dataDir + "angles" + dn + ".txt");
thrust  = readmatrix(dataDir + "thrust" + dn + ".txt");
time    = readmatrix(dataDir + "time" + dn + ".txt");
kept    = readmatrix(dataDir + "kept" + dn + ".txt");
consts  = readmatrix(dataDir + "consts" + dn + ".txt");
sunangs = readmatrix(dataDir + "sunangles" + dn + ".txt");

% Initial keplerian elements
kep0    = [24363.99406945558;
           0.7306;
           0.4974187495004638;
           0.0;
           0.0;
           0.0];

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
tmax   = max(thrust);
sf     = zeros(size(cart, 1), 1);
dir    = nan(size(cart, 1), 3);
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

    % Compute thrust direction
    a = angles(i,1);
    b = angles(i,2);
    sf(i) = thrust(i) / tmax;
    dir(i,:) = [cos(b)*sin(a), cos(b)*cos(a), sin(b)];
end

% Construct output file for video creation
output = [time, cart, dir, sf];
writematrix(output, "data/gto_to_geo_trajectory.csv");

% 3D Trajectory
figure()
plot3(cart_t(:,1),cart_t(:,2),cart_t(:,3),"Color","#FF0000")
hold on
plot3(cart_c(:,1),cart_c(:,2),cart_c(:,3),"b")
plot3(cart_e(:,1),cart_e(:,2),cart_e(:,3),"g")
plotKeplerianOrbit(kept, 360, "3D")
xlabel("x")
ylabel("y")

axis equal
grid on

% X-Y Proj. of Trajectory
figure()
plot(cart_t(:,1),cart_t(:,2),"Color","#E56A54")
hold on
plot(cart_e(:,1),cart_e(:,2),"Color","#005bbb")
plotKeplerianOrbit(kep0, 360, "XY")
plotKeplerianOrbit(kept, 360, "XY")
xlabel("x, km")
ylabel("y, km")
%leg111 = legend("Thrusting", "Eclipsed","location","northwestoutside");

axis equal
padding = 5e3;
%xlim([-15e4,max([max(cart_t(:,1)),max(cart_e(:,1))])] + padding)
xlim([min([min(cart_t(:,1)),min(cart_e(:,1))]) - padding, ...
      max([max(cart_t(:,1)),max(cart_e(:,1))]) + padding])
ylim([min([min(cart_t(:,2)),min(cart_e(:,2))]) - padding, ...
      max([max(cart_t(:,2)),max(cart_e(:,2))]) + padding])
grid on

% Publication plotting settings
set(gca, "fontname", "Times New Roman", "fontsize", 10)
%set(leg111, "fontname", "Times New Roman", "fontsize", 8)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,3,3])
set(gcf, "PaperPositionMode","Manual")

% Print figure
print(figFolder + "traj_xy_" + dn + ".png", "-dpng", "-r300")

figure()
tiledlayout(3,1)

% X-Y Plot
nexttile([2,1])
plot(cart_t(:,1),cart_t(:,2),"r")
hold on
plot(cart_e(:,1),cart_e(:,2),"g")
plot(cart_c(:,1),cart_c(:,2),"b")
plotKeplerianOrbit(kep0, 360, "XY")
plotKeplerianOrbit(kept, 360, "XY")

% Setup axes
xlabel("x, km", "Interpreter", "Latex")
ylabel("y, km", "Interpreter", "Latex")
axis equal
grid on

% Setup legend
leg1 = legend("Thrusting","Eclipsed","Coasting","Location","northwest");

% Publication plotting settings
set(gca, "fontname", "Times New Roman", "fontsize", 10)

% Y-Z Plot
nexttile
plot(cart_t(:,2),cart_t(:,3),"r")
hold on
plot(cart_e(:,2),cart_e(:,3),"g")
plot(cart_c(:,2),cart_c(:,3),"b")
plotKeplerianOrbit(kep0, 360, "YZ")
plotKeplerianOrbit(kept, 360, "YZ")

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
    else
        i = j;
    end
    nexttile
    plot(time, kep(:,i))
    hold on
    if j == 1
        yline(42165.0,"k")
        ylabel("a, km", "Interpreter", "Latex")
        ylim([0.0, 6.7e4])
    elseif j == 2
        yline(0.0,"k")
        ylabel("e", "Interpreter", "Latex")
        ylim([-0.1, 0.85])
        yticks([0,0.3,0.6])
    elseif j == 3
        yline(0.0,"k")
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

for i = 1:length(sunangs)
    if sunangs(i) > 90
        sunangs(i) = 180 - sunangs(i);
    end
end

figure()
scatter(time, sunangs,0.5,".", "MarkerEdgeColor", "#005bbb", "MarkerFaceColor", "#005bbb")
hold on
%if dn == "cdp"
%    yline([90 - 40, 90 + 40],"k")
%else
%    yline(90 - 40,"k")
%end
yline(50, "k")
xlabel("time, days", "Interpreter", "Latex")
%ylabel("$\psi$, deg", "Interpreter", "Latex")
ylabel("$\phi$, deg", "Interpreter", "Latex")
grid on

% Publication plotting settings
set(gca, "fontname", "Times New Roman", "fontsize", 10)
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,3.0,3.0])
set(gcf, "PaperPositionMode","Manual")

% Print figure
%print(figFolder + "sunAngles_" + dn + ".eps", "-depsc", "-r300")
%print(figFolder + "sunAngles_" + dn + ".eps", "-depsc", "-r300")
print(figFolder + "sunAngles_" + dn + ".png", "-dpng", "-r300")

% FOR FOUR-SQUARE ONLY
set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,5.0,4.0])
print(figFolder + "fourSquareAngles.pdf", "-dpdf","-r300")

% % Zoomed in Sun Angle Plot
% figure()
% scatter(time,sunangs,".")
% hold on
% if dn == "cdp"
%     yline([90 - 40, 90 + 40],"k")
% else
%     yline(90 - 40,"k")
% end
% xlim([0.2543,0.6447])
% xlabel("time, days", "Interpreter", "Latex")
% ylabel("$\psi$, deg", "Interpreter", "Latex")
% grid on
% 
% % Publication plotting settings
% set(gca, "fontname", "Times New Roman", "fontsize", 10)
% set(gcf, "PaperUnits","inches","PaperPosition",[0.25,0.25,5.0,1.5])
% set(gcf, "PaperPositionMode","Manual")
% 
% % Print figure
% print(figFolder + "sunAnglesZoomed_" + dn + ".eps", "-depsc", "-r300")

