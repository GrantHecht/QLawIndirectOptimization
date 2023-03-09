clear; close all; clc

% Define symbolic variables
syms a1 a2 a3 real
a = [a1; a2; a3];

syms u1 u2 u3 real
u = [u1; u2; u3];

syms lam1 lam2 real
syms gamma amax real positive

% Define required matricies
C = [1, 0, 0; 0, 1, 0; 0, 0, 0];
d = [0; 0; 1.0 / gamma];

% ===== Case 1: -a lies outside of or on cone
u1 = amax*[-a1; -a2; gamma*sqrt(a1^2 + a2^2)] / ...
        sqrt((1+gamma^2)*a1^2 + (1+gamma^2)*a2^2);

% Primal feasibility
pf11    = simplify(norm(u1) - amax)
pf12    = simplify(norm(C*u1) - d'*u1)

% Stationary condition and dual feasibility
gL1     = a + lam1*u1/norm(u1) + lam2*(C*u1/norm(C*u1) - d);
lams1   = solve(gL1 == zeros(3,1), [lam1,lam2]);
fprintf("Dual Feasibility Requirements:\n")
lam11   = simplify(lams1.lam1)
lam12   = simplify(lams1.lam2)

% ===== Case 2: -a lies inside of cone
u2 = -amax*a/norm(a);

% Primal feasibility
pf21    = simplify(norm(u2) - amax)
pf22    = simplify(norm(C*u2) - d'*u2)

% Stationary condition and dual feasibility
gL2     = a + lam1*u2/norm(u2);
lam21   = simplify(solve(gL2 == zeros(3,1), lam1))

