clear; close all; clc

% Define symbolic variables
syms a1 a2 a3 real
a = [a1; a2; a3];

syms u1 u2 u3 real
u = [u1; u2; u3];

syms lam1 lam2 real
syms gamma amax Tmax m real positive

% Define required matricies
C = [1, 0, 0; 0, 1, 0; 0, 0, 0];
d = [0; 0; 1.0 / gamma];

% % ===== Case 1: -a lies outside of or on cone
% u1 = amax*[-a1; -a2; gamma*sqrt(a1^2 + a2^2)] / ...
%         sqrt((1+gamma^2)*a1^2 + (1+gamma^2)*a2^2);
% u1t = amax*[-a1; -a2; gamma*sqrt(a1^2 + a2^2)] / ...
%         sqrt((1+gamma^2)*(a1^2 + a2^2));
% 
% % Primal feasibility
% pf11    = simplify(norm(u1) - amax)
% pf12    = simplify(norm(C*u1) - d'*u1)
% 
% % Stationary condition and dual feasibility
% gL1     = a + lam1*u1/norm(u1) + lam2*(C*u1/norm(C*u1) - d);
% lams1   = solve(gL1 == zeros(3,1), [lam1,lam2]);
% fprintf("Dual Feasibility Requirements:\n")
% lam11   = simplify(lams1.lam1)
% lam12   = simplify(lams1.lam2)
% 
% % ===== Case 2: -a lies inside of cone
% u2 = -amax*a/norm(a);
% 
% % Primal feasibility
% pf21    = simplify(norm(u2) - amax)
% pf22    = simplify(norm(C*u2) - d'*u2)
% 
% % Stationary condition and dual feasibility
% gL2     = a + lam1*u2/norm(u2);
% lam21   = simplify(solve(gL2 == zeros(3,1), lam1))

% Case 4:
% Stationary condition
sc4     = a + lam1*(m/Tmax)*u + lam2*((C*u)/(d'*u) - d);

% Substitute ucz = gamma*sqrt(ucx^2 + ucy^2) into stationary condition
sc4_ucz = subs(sc4, u3, gamma*sqrt(u1^2 + u2^2));

% Solve simplified stationary condition to get all a on lhs
cond1 = isolate(simplify(sc4_ucz(1)) == 0, u1)
cond2 = isolate(simplify(sc4_ucz(2)) == 0, u2)
cond3 = simplify(subs(sc4_ucz(3)) == 0)

% b86     = (lam2*(gamma^2 + 1) - a3*gamma)/(gamma^2*sqrt(u1^2 + u2^2));
% b1      =  (m/Tmax)*sqrt((gamma^2 + 1)*(a1^2 + a2^2));
% b2      = -(m/Tmax)*sqrt((gamma^2 + 1)*(a1^2 + a2^2));
% 
% lam21   = simplify(solve(b86 == b1, lam2))
% lam22   = simplify(solve(b86 == b2, lam2))
