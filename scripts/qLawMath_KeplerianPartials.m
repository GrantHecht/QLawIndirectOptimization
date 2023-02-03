%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script calculates the derivatives of the Lyapunov function of
% Petropoulos's control law with respect to semi-major axis,
% eccentricity, inclination, right ascension of the ascending
% node, and argument of periapsis symbolically.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This script uses the same variable names as Nobel Hatten's script from his
%master's thesis
%Nobel's thesis does not include "Optimal" control formulation and C code
%generation, these features were added by Donald Ellison Feb. 6th 2014 
% Modifications added by Grant Hecht starting Jan. 27th 2023
%(see note below)

% Declare symbolic variables
syms Q dQdt real    % Lyapunov function and its time derivative
syms mu real        % gravitational parameter of central body
syms f real         % magnitude of control acceleration
syms fr ftheta fh real  % components of control acceleration in R-theta-H system
syms alpha beta real
syms h p r real % osculating values of: angular momentum magnitude, semi-latus rectum, and position vector magnitude

syms Wp P rper rpermin real % penalty function variables

syms W Wsma We Winc Wape Wran real  %COE selection weights

syms sma e inc ape ran tru real         % osculating orbital elements
syms sma_t e_t inc_t ape_t ran_t real   % target orbital elements
syms Qsma Qe Qinc Qape Qran real        % terms of Lyapunov function summation
syms smadotxx edotxx incdotxx randotxx apedotxxi apedotxxo real
syms apedotxx cosvxxo rxxo distape distran m_petro n_petro real
syms r_petro b_petro k_petro real

% ^^^ maximum rates of change of each element over
% osculating orbit and other variables used in their calculation
syms dQdcoe real% 5x5 array holding derivatives of each term of

% Lyapunov function with respect to each orbital element
% put h, p, r in terms of orbital elements
p = sma*(1-e^2);
h = sqrt(mu*p);
r = p/(1+e*cos(tru));
rper = sma*(1-e);

% maximum values of rates of change of orbital elements over
% control acceleration angles and true anomaly

% semi-major axis
smadotxx = 2*f*sqrt((sma^3*(1+e))/(mu*(1-e)));

% eccentricity
edotxx = (2*p*f)/h;

% inclination
incdotxx = (p*f)/(h*(sqrt(1-e^2*(sin(ape))^2)-e*abs(cos(ape))));

% right ascension of the ascending node
randotxx = (p*f)/(h*sin(inc)*(sqrt(1-e^2*(cos(ape))^2)- ...
e*abs(sin(ape))));

% argument of periapsis: in-plane
cosvxxo = (((1-e^2)/(2*e^3)) + ...
sqrt((1/4)*((1-e^2)/e^3)^2+(1/27)))^(1/3) - ...
(((-(1-e^2)/(2*e^3))) + ...
sqrt((1/4)*((1-e^2)/e^3)^2+(1/27)))^(1/3) - (1/e);
rxxo = p/(1+e*cosvxxo);
apedotxxi = (f/(e*h))*sqrt(p^2*cosvxxo^2 + (p+rxxo)^2* ...
(1-cosvxxo^2));

% argument of periapsis: out-of-plane
apedotxxo = randotxx*abs(cos(inc));

% Minimum periapse penalty function
P = exp(k_petro*(1 - rper/rpermin));

% build Lyapunov function term by term

% semi-major axis term also has scaling function:
% Grant Hecht - Removed absolute value from equation
s_sma = (1+((sma-sma_t)/(m_petro*sma_t))^n_petro)^(1/r_petro);

% semi-major axis term
Qsma = s_sma*(((sma - sma_t)/smadotxx)^2);

% eccentricity term
Qe = ((e - e_t)/edotxx)^2;

% inclination term
Qinc = ((inc - inc_t)/incdotxx)^2;

% right ascension of the ascending node term
distran = acos(cos(ran - ran_t));
Qran = (distran/randotxx)^2;

% argument of periapsis term
distape = acos(cos(ape - ape_t));
apedotxx = (apedotxxi+b_petro*apedotxxo)/(1+b_petro);
Qape = (distape/apedotxx)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This point on, additions by Donald Ellison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COE selection weights
W = [Wsma We Winc Wape Wran];

% Build Q
Q = (1 + Wp*P)*(Wsma*Qsma + We*Qe + Winc*Qinc + Wape*Qape + Wran*Qran);

% partial differentials of full Q with penalty term
dQfulldcoe(1,1) = diff(Q, sma);
dQfulldcoe(1,2) = diff(Q, e);
dQfulldcoe(1,3) = diff(Q, inc);
dQfulldcoe(1,4) = diff(Q, ape);
dQfulldcoe(1,5) = diff(Q, ran);

% Generate code for partials
ccode(dQfulldcoe(1,1), 'file', './data/dQdsma.txt');
ccode(dQfulldcoe(1,2), 'file', './data/dQde.txt');
ccode(dQfulldcoe(1,3), 'file', './data/dQdinc.txt');
ccode(dQfulldcoe(1,4), 'file', './data/dQdape.txt');
ccode(dQfulldcoe(1,5), 'file', './data/dQdran.txt');
ccode(dQfulldcoe, 'file', './data/Qpartials.txt');



