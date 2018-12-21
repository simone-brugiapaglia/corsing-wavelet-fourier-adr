function A = stiffmatWF_ADR2D_single(L,etas,betas,rhos,tau_grid,opt)
% A = STIFFMATWF_ADR2D_SINGLE(L,etas,betas,rhos,tau_grid,opt) computes the 
% stiffness matrix of the Petrov-Galerkin discretization of an
% advection-diffusion-reaction equation with separable coefficients
% associated with a sigle-scale trial basis and a selection of Fourier test
% functions. 
%
% L     = max level (mesh size h = 2^(-L))
% etas  = {eta1,  eta2}
% betas = {beta1, beta2}
% rhos  = {rho1,  rho2}
% tau_grid  = tau_grid(1,:) x tau_grid(2,:) are the test indices
% opt   = options struct as in stiffmat_ADR1D_single
%
% A = stiffness matrix

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

OO = @(x) 0*x;

if (~isconst(etas{1}) || ~isconst(etas{2}) || ~isconst(betas{1,1}) || ...
        ~isconst(betas{1,2}) || ~isconst(betas{2,1}) || ...
        ~isconst(betas{2,2}) || ~isconst(rhos{1}) || ~isconst(rhos{2}) )
    error('2D implementation tensted only for constant coefficients')
end

%% Compute 1D stiffness matrices
K_eta  = cell(1,2);  % {K_eta_1, K_eta_2}
M_eta  = cell(1,2);  % {M_eta_1, M_eta_2}
M_rho  = cell(1,2);  % {M_rho_1, M_rho_2}
T_beta = cell(1,2); % {T_beta_11, T_beta_22}
M_beta = cell(1,2); % {M_beta_12, M_beta_21}
% We use the decomposition 
% B = K_eta_1 x M_eta_2 + M_eta_1 x K_eta_2 
%   + T_beta_11 x M_beta_12 + M_beta_12 x T_beta_22
%   + M_rho_1 x M_rho_2,
% where are discretization of 1D problems defined as follows:
% (K_f)_ij = (f psi_j', xi_i') (diffusion)
% (T_f)_ij = (f psi_j', xi_i)  (advection)
% (M_f)_ij = (f psi_j, xi_i)   (reaction)
for i_dim = 1:2
    K_eta{i_dim}  = stiffmatWF_ADR1D_single(L,etas{i_dim},OO,OO,tau_grid(i_dim,:),opt);
    M_eta{i_dim}  = stiffmatWF_ADR1D_single(L,OO,OO,etas{i_dim},tau_grid(i_dim,:),opt);
    M_rho{i_dim}  = stiffmatWF_ADR1D_single(L,OO,OO,rhos{i_dim},tau_grid(i_dim,:),opt);
    T_beta{i_dim} = stiffmatWF_ADR1D_single(L,OO,betas{i_dim,i_dim},OO,tau_grid(i_dim,:),opt);
    i_next = mod((i_dim + 1) - 1, 2) + 1; % if i_dim = 1, i_next = 2; if i_dim = 2, i_next = 1
    M_beta{i_dim} = stiffmatWF_ADR1D_single(L,OO,OO,betas{i_dim,i_next},tau_grid(i_dim,:),opt);    
end

%% Assemble the 2D stiffness matrix
A = kron(K_eta{1},M_eta{2}) + kron(M_eta{1},K_eta{2}) + ...
    kron(T_beta{1},M_beta{1}) + kron(M_beta{2},T_beta{2}) + ...
    kron(M_rho{1},M_rho{2});


end

