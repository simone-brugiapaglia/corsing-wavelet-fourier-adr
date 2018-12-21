function [f1,f2] = forcingterm2D(u1,u2,etas,betas,rhos)
% F = FORCINGTERM2D(U,ETA,BETA,RHO) computes the forcing term of a
% two-dimensional advection-diffusion-reaction equation with separable 
% coefficients defined by ETAS, BETAS, and RHOS by symbolic computations 
% assuming that the exact solution is sum of separable terms stored in
% U1 and U2.
%
% EAT  = etas{1} x etas{2};
% BETA = [betas{1,1} x betas{1,2};
%         betas{2,1} x betas{2,2}];
% RHO = rhos{1} x rhos{2};
%
%        K                            5K
%       ---                           ---
% U   = >   U1{i} x U2{i}       F   = >   F1{i} x F2{i}
%       ---                           ---
%       i=1                           i=1
%
%

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

%% Symbolic computations
syms x;
f1 = cell(1,length(u1)*5);
f2 = cell(1,length(u1)*5);
for i_deco = 1:length(u1)
    f1{5 * (i_deco-1) + 1} = matlabFunction( - diff(etas{1}(x) * diff(u1{i_deco}(x))) );
    f2{5 * (i_deco-1) + 1} = matlabFunction( etas{2}(x) * u2{i_deco}(x) );
    
    f1{5 * (i_deco-1) + 2} = matlabFunction( etas{1}(x) * u1{i_deco}(x) );
    f2{5 * (i_deco-1) + 2} = matlabFunction( - diff(etas{2}(x) * diff(u2{i_deco}(x))) );
    
    f1{5 * (i_deco-1) + 3} = matlabFunction( betas{1,1}(x) * diff(u1{i_deco}(x)) );
    f2{5 * (i_deco-1) + 3} = matlabFunction( betas{1,2}(x) * u2{i_deco}(x) );
    
    f1{5 * (i_deco-1) + 4} = matlabFunction( betas{2,1}(x) * u1{i_deco}(x) );
    f2{5 * (i_deco-1) + 4} = matlabFunction( betas{2,2}(x) * diff(u2{i_deco}(x)) );
    
    f1{5 * (i_deco-1) + 5} = matlabFunction( rhos{1}(x) * u1{i_deco}(x) );
    f2{5 * (i_deco-1) + 5} = matlabFunction( rhos{2}(x) * u2{i_deco}(x) );
end
clear x;

%% Conversion of constant function handles
% convert handles of the form " @() c " to " @(x) 0*x + c "
for i = 1:length(f1)
    if nargin(f1{i}) == 0
        f1{i} = @(x) 0*x + f1{i}();
    end
    if nargin(f2{i}) == 0
        f2{i} = @(x) 0*x + f2{i}();
    end
end
