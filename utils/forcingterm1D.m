function f = forcingterm1D(u,eta,beta,rho)
% F = FORCINGTERM1D(U,ETA,BETA,RHO) computes the forcing term of and
% advection-diffusion-reaction equation with coefficients ETA, BETA, and
% RHO by symbolic computations assuming that the exact solution is U.


% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)
syms x;
% symbolic differentiation
fSym = - diff(eta(x) * diff(u(x))) + beta(x) * diff(u(x)) + rho(x) * u(x); 
f = matlabFunction(fSym); 
