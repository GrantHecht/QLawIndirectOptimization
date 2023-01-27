clear; close all; clc

% Define variables
syms p f g h k L real % Modified Eq. Elements
syms a e inc ran aop nu real % Kep. Elemenents

% Make assumptions
assume(e > 0 & e < 1)
assume(inc > 0 & inc < pi)
assume(a > 0)
assume(p > 0)
assume(-1 < f & f < 1)
assume(-1 < g & g < 1)

eqs = [p == a*(1 - e^2);
       f == e*cos(aop + ran);
       g == e*sin(aop + ran);
       h == tan(inc/2)*cos(ran);
       k == tan(inc/2)*sin(ran);
       L == ran + aop + nu];


seqs = solve(eqs, [a,e,inc,ran,aop,nu]);


