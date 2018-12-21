function A = stiffmatWF_ADR1D_multi(l0,L,eta,beta,rho,tau,opt)
% BUILD_STIFFMAT builds the stiffness matrix associated with the
% Petrov-Galerkin discretization of an advection-diffusion-reaction
% problem on the torus T = R/Z. 
%
% In particular, the bilinear form associate with the weak formulation of 
% the problem is
%
% a(u,v) = ( eta Du, Dv ) + ( beta Du, v ) + ( rho u, v );
%
% where ( u,v ) = int_T u conj(v).
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
%  opt   .forcequad  = 1 to force numerical quadrature 0 otherwise
%        .L_add  = additional levels to improve quadrature accuracy
%        .trial  = 'bior2.2'  (default) Biortogonal B-splines d = dtilde = 2
%                = 'hier'     Hierarchical Multisale Basis
%               
%
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
% EPFL, 2016 -- SFU, 2018

if nargin == 6
    opt.forcequad = 0;
    opt.L_add = 0;
    opt.trial = 'bior';
else
    if ~isfield(opt,'forcequad')
        opt.forcequad = 0;
    end
    if ~isfield(opt,'L_add')
        opt.L_add = 0;
    end
    if ~isfield(opt,'trial')
        opt.trial = 'bior';
    end
    
end

N = 2^L;

% Wavelet filters
[primal,dual] = filterbank(opt.trial);

L_add = opt.L_add;
A_single = stiffmatWF_ADR1D_single(L+L_add,eta,beta,rho,tau,opt);

% Single-scale to multi-scale conversion
A_aug = wavedeco1D(A_single,l0,L+L_add,primal);

% Select only the levels up to L
A1 = A_aug(:,1:N);

%% Normalization wrt H1-norm   
DXi  = diag(1./sqrt(1 + (2*pi*tau).^2));
ll = max(l0,floor(log2(0:N-1)));
DPsi = diag(2.^(-ll));

A = DXi * A1 * DPsi; 
