function A = stiffmatWF_ADR1D_single(L,eta,beta,rho,tau,opt)
% BUILD_STIFFMAT builds the stiffness matrix associated with the
% Petrov-Galerkin discretization of an advection-diffusion-reaction
% problem on the torus T = R/Z.
%
% In particular, the bilinear form associate with the weak formulation of
% the problem is
%
% a(u,v) = < eta Du, Dv > + < beta Du, v > + < rho u, v >;
%
% where < u,v > = int_{T} u conj(v).
%
% The trial functions are the biorthogonal B-spline wavelets of parameters
% d = dtilde = 2 (see [1]). The test functions are a selection of the Fourier basis
% functions {exp(2i pi q x) : q in Z}, normalized wrt the H^1(T)-norm.
%
%
%
% INPUT:
%  l0    coarsest level
%  L     finest level
%  eta   diffusion term (function handle)
%  beta  advection term (function handle)
%  rho   reaction term (function handle)
%  tau   test indices
%  opt   .forcequad  == 1 force the usage of quadrature formulas
%
% OUTPUT:
%  A     stiffness matrix with entries a(psi_{l,k},xi_{tau(i)})
%
%
% REFERENCES:
%
% [1] K. Urban. Wavelet Methods for Elliptic Partial Differential
% Equations. Oxford University Press, 2009.

% Simone Brugiapaglia, Fabio Nobile
% EPFL, 2016

if nargin == 5
    opt.forcequad = 0;
end


% define relevant parameters
N = 2^L;
m = length(tau);
h = 2^(-L);

if (isconst(eta) && isconst(beta) && isconst(rho) && ~(opt.forcequad))
    %% CONSTANT COEFFICIENTS CASE
    
    % Initialization
    A = zeros(m,N);
    
    % Build A row by row
    for i = 1:m
        % Explicit formula for the diffusive term
        q = tau(i);
        kk = 0:N-1;
        A_row_di = 4* 2^(3*L/2) .* exp(-2i*pi*q*kk*h) .* ...
            sin(pi* q * ones(1,N) * h).^2;
        if tau(i) == 0
            A_row_ad = zeros(1,N);
            A_row_re = 2^(-L/2) * ones(1,N);
        else
            % exploit integration by parts and the fact that 
            % xi_q' = (2i pi q) xi_q to compute the advective and reactive
            % terms
            A_row_ad = -1/(2i*pi*q)  * A_row_di;
            A_row_re = 1/(2i*pi*q) * A_row_ad;
        end
        
        A(i,:) = (eta(0)*A_row_di + beta(0)*A_row_ad + rho(0)*A_row_re);
    end
    

    
else

    
    %% NONCONSTANT COEFFICIENTS CASE
    
    
    % grid on the torus T
    xh = 0:h:1-h;
    
    % build the mass matrix M_kl = (phi_{L,k},phi_{L,l})
    v = zeros(1,N);
    v(1) = 4;
    v(2) = 1;
    v(end) = 1;
    M = sparse(1/6 * toeplitz(v));
    
    % initialization
    A_single_ad = zeros(m,N);
    A_single_r  = zeros(m,N);
    
    % In this for loop we compute the Petrov-Galerkin discretization associated
    % with the single scale trial functions, i.e., the family of hat functions
    % associated with the grid of step h on the torus.
    for i = 1:m
        % define current test function xi and its derivative
        q = tau(i);
        xi  = @(x) exp(2i*pi*q*x);
        Dxi = @(x) (2i*pi*q) * exp(2i*pi*q*x);
        
        % advection-diffusion term
        ad_fun = @(x) eta(x) .* conj(Dxi(x)) + beta(x) .* conj(xi(x));
        
        % we approximate the term
        %
        %                  int_T ad_fun Dphi_{L,k},
        %
        % i.e., the integral of ad_fun times the derivative of the k-th hat
        % function at the finest level L, employing the composite trapezoidal
        % qudrature rule of step h.
        % We notice that Dphi_{L,k} is a suitably rescaled step function,
        % supported on the interval [h(k-1), h(k+1)], with a jump on the
        % midpoint.
        
        % % Trapezoidal quadrature rule
        xh_l = circshift(xh,[0,1]);
        xh_r = circshift(xh,[0,-1]);
        % In order to derive the following formula, recall that Dphi_{L,k} is a
        % step function rescaled by a factor 2^(3L/2).
        A_single_ad(i,:) = 2^(L/2)/2 * (ad_fun(xh_l) - ad_fun(xh_r));
        
        
        % the reaction term, i.e.,
        %
        %                ctilde(k) = int_T r_fun phi_{L,k},
        %
        % is approximated computing the coefficients of r_fun with respect to
        % the biorthogonal basis {phitilde_{L,k}} of {phi_{L,k}}. Indeed, we
        % have
        %
        %                r_fun = sum ctilde(k) phitilde_{L,k}
        %
        % ctilde is computed in two steps:
        %
        % 1) compute the coefficients {c(k)} of r_fun with respect to the primal
        % basis by evaluation at the nodal points and rescaling by a factor
        % 2^(-L/2). This is due to the fact that {2^(-L/2) phi_{L,k}} is
        % lagrangian with respect to the uniform grid.
        %
        % 2) map the primal coefficients {c(k)} into the dual coefficients
        % {ctilde(k)} thorugh the mass matrix M, defined as
        %
        %               M(k,j) = int_T  phi_{L,k} phi_{L,j}.
        r_fun = @(x) rho(x).*conj(xi(x));
        rh    = r_fun(xh);
        A_single_r(i,:)  = (rh/(2^(L/2))) * M;
    end

    A = A_single_ad + A_single_r; % (Advection-diffusion) + (reaction)
    
end
