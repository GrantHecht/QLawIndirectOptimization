clear; close all; clc

% Read in data
lams    = readmatrix("data/TaheriGTO2GEO_costates_new.txt");
meefs   = readmatrix("data/TaheriGTO2GEO_finalStates_new.txt");
tfs     = readmatrix("data/TaheriGTO2GEO_tofs_new.txt");
percs   = readmatrix("data/TaheriGTO2GEO_percs_new.txt");
ftols   = readmatrix("data/TaheriGTO2GEO_ftols_new.txt");

% Plot
tiledlayout(7,1)

for i = 1:6
    nexttile 
    plot(lams(:,i))
    if i == 1
        ylabel("$\lambda_p$", "Interpreter","latex")
    elseif i == 2
        ylabel("$\lambda_f$", "Interpreter","latex")
    elseif i == 3
        ylabel("$\lambda_g$", "Interpreter","latex")
    elseif i == 4
        ylabel("$\lambda_h$", "Interpreter","latex")
    elseif i == 5
        ylabel("$\lambda_k$", "Interpreter","latex")
    else
        ylabel("$\lambda_L$", "Interpreter","latex")
    end
end
nexttile 
semilogy(ftols)