function [f1,f2,f3] = forcingterm3D(u1,u2,u3,etas,betas,rhos)
% [f1,f2,f3] = FORCINGTERM3D(u1,u2,u3,etas,betas,rhos) computes the forcing
% term of the advection-diffusion-reaction equation given separable 
% coefficients eta, beta, rho, and the exact solution in low-rank form.
%
% Input:
% u1, u2, u3, etas, betas, rhos = cell array such that
%
% u    = sum_i u1{i} x u2{i} x u3{i} 
% eta  = etas{1} x etas{2} x etas{3}
% beta = [betas{1,1} x betas{1,2} x betas{1,3};
%         betas{2,1} x betas{2,2} x betas{2,3};
%         betas{3,1} x betas{3,2} x betas{3,3};
% rho  = rhos{1} x rhos{2} x rhos{3}
%
% Output:
% f1, f2, f3 = cell arrays such that
%
% f    = sum_i f1{i} x f2{i} x f3{i}; 

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

if length(u1) ~= length(u2) || length(u1) ~= length(u3)
    error('u1, u2, and u3 must have the same length.')
end

%% Symbolic computation of the forcing term
syms x;
f1 = cell(1,length(u1)*7);
f2 = cell(1,length(u1)*7);
f3 = cell(1,length(u1)*7);
for i_deco = 1:length(u1)
    % Diffusion 1: -D(mu1 D(u1)) (mu2 u2) (mu3 u3)
    f1{7 * (i_deco-1) + 1} = matlabFunction( - diff(etas{1}(x) * diff(u1{i_deco}(x))) );
    f2{7 * (i_deco-1) + 1} = matlabFunction(        etas{2}(x) *      u2{i_deco}(x) );
    f3{7 * (i_deco-1) + 1} = matlabFunction(        etas{3}(x) *      u3{i_deco}(x) );
    
    % Diffusion 2: (mu1 u1) [-D(mu2 D(u2))] (mu3 u3)
    f1{7 * (i_deco-1) + 2} = matlabFunction(        etas{1}(x) *      u1{i_deco}(x) );
    f2{7 * (i_deco-1) + 2} = matlabFunction( - diff(etas{2}(x) * diff(u2{i_deco}(x))) );
    f3{7 * (i_deco-1) + 2} = matlabFunction(        etas{3}(x) *      u3{i_deco}(x) );
    
    % Diffusion 3: (mu1 u1) (mu2 u2) [-D(mu3 D(u3))]
    f1{7 * (i_deco-1) + 3} = matlabFunction(        etas{1}(x) *      u1{i_deco}(x) );
    f2{7 * (i_deco-1) + 3} = matlabFunction(        etas{2}(x) *      u2{i_deco}(x) );
    f3{7 * (i_deco-1) + 3} = matlabFunction( - diff(etas{3}(x) * diff(u3{i_deco}(x))) );
    
    % Advection 1: (beta11 D(u1)) (beta12 u2) (beta13 u3)    
    f1{7 * (i_deco-1) + 4} = matlabFunction( betas{1,1}(x) * diff(u1{i_deco}(x)) );
    f2{7 * (i_deco-1) + 4} = matlabFunction( betas{1,2}(x) *      u2{i_deco}(x) );
    f3{7 * (i_deco-1) + 4} = matlabFunction( betas{1,3}(x) *      u3{i_deco}(x) );
    
    % Advection 2: (beta21 u1) (beta22 D(u2)) (beta23 u3)
    f1{7 * (i_deco-1) + 5} = matlabFunction( betas{2,1}(x) *      u1{i_deco}(x) );
    f2{7 * (i_deco-1) + 5} = matlabFunction( betas{2,2}(x) * diff(u2{i_deco}(x)) );
    f3{7 * (i_deco-1) + 5} = matlabFunction( betas{2,3}(x) *      u3{i_deco}(x) );
    
    % Advection 3: (beta31 u1) (beta32 u2) (beta33 D(u3))
    f1{7 * (i_deco-1) + 6} = matlabFunction( betas{3,1}(x) *      u1{i_deco}(x) );
    f2{7 * (i_deco-1) + 6} = matlabFunction( betas{3,2}(x) *      u2{i_deco}(x) );
    f3{7 * (i_deco-1) + 6} = matlabFunction( betas{3,3}(x) * diff(u3{i_deco}(x)) );
    
    % (rho1 u1) (rho2 u2) (rho3 u3)
    f1{7 * (i_deco-1) + 7} = matlabFunction( rhos{1}(x) * u1{i_deco}(x) );
    f2{7 * (i_deco-1) + 7} = matlabFunction( rhos{2}(x) * u2{i_deco}(x) );
    f3{7 * (i_deco-1) + 7} = matlabFunction( rhos{3}(x) * u3{i_deco}(x) );
end
clear x;


% convert handles of the form " @() c " to " @(x) 0*x + c "
for i = 1:length(f1)
    if nargin(f1{i}) == 0
        f1{i} = @(x) 0*x + f1{i}();
    end
    if nargin(f2{i}) == 0
        f2{i} = @(x) 0*x + f2{i}();
    end
    if nargin(f3{i}) == 0
        f3{i} = @(x) 0*x + f3{i}();
    end
end