function B = stiffmatWF_ADR3D_single(L,etas,betas,rhos,tau_grid,opt)
% A = STIFFMATWF_ADR3D_SINGLE(L,etas,betas,rhos,tau_grid,opt) computes the 
% stiffness matrix of the Petrov-Galerkin discretization of an
% advection-diffusion-reaction equation with separable coefficients
% associated with a sigle-scale trial basis and a selection of Fourier test
% functions. 
%
% L     = max level (mesh size h = 2^(-L))
% etas  = {eta1,  eta2, eta3}
% betas{i,j} = beta_ij
% rhos  = {rho1,  rho2, rho3}
% tau_grid  = tau_grid(1,:) x tau_grid(2,:) x tau_grid(3,:) are the test indices
% opt as in stiffmatWF_ADR1D
%
% A = stiffness matrix
%
% a(u,v) = (eta * grad(u), grad(v)) + (beta * grad(u) , v) + (rho * u , v),
%
% where the coefficients are separable 
%  eta(x,y)  = eta_1(x) eta_2(y) eta_3(z),
%  beta(x,y) = [beta_11(x), beta_12(y), beta_13(z) ; 
%               beta_21(x), beta_22(y), beta_23(z) ;
%               beta_31(x), beta_32(y), beta_33(z) ],
%  rho(x,y)  = rho_1(x) rho_2(y) rho_3(z).

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

if (~isconst(etas{1}) || ~isconst(etas{2}) || ~isconst(etas{3}) || ...
        ~isconst(betas{1,1}) || ~isconst(betas{1,2}) || ~isconst(betas{1,3}) || ...
        ~isconst(betas{2,1}) || ~isconst(betas{2,2}) || ~isconst(betas{2,3}) || ...
        ~isconst(betas{3,1}) || ~isconst(betas{3,2}) || ~isconst(betas{3,3}) || ...
        ~isconst(rhos{1}) || ~isconst(rhos{2}) || ~isconst(rhos{3}) )
    warning('3D stiffness matrix''s accuracy may be low for nonconstant coefficients')
end

OO = @(x) 0*x;



% assemble 1D-stiffness matrices
K_eta = cell(1,3);  % {K_eta_1, K_eta_2, K_eta_3}
M_eta = cell(1,3);  % {M_eta_1, M_eta_2, M_eta_3}
M_rho = cell(1,3);  % {M_rho_1, M_rho_2, M_rho_3}
T_beta = cell(1,3); % {T_beta_11, T_beta_22, T_beta_33}
M_beta = cell(3,2); % {M_beta_12, M_beta_13; 
                    %  M_beta_21, M_beta_23; 
                    %  M_beta_31, M_beta_32}
% We use the decomposition 
% B = K_eta_1 x M_eta_2 X M_eta_3 + M_eta_1 x K_eta_2 x M_eta_3 + M_eta_1 x M_eta_2 x K_eta_3 
%   + T_beta_11 x M_beta_12 x M_beta_13 + M_beta_12 x T_beta_22 x M_beta_23 + M_beta_31 x M_beta_32 x T_beta_33   
%   + M_rho_1 x M_rho_2 x M_rho_3 ,
% where are discretization of 1D problems defined as follows:
% (K_f)_ij = (f psi_j', xi_i') (diffusion)
% (T_f)_ij = (f psi_j', xi_i)  (advection)
% (M_f)_ij = (f psi_j, xi_i)   (reaction)
for i_dim = 1:3
    K_eta{i_dim}  = stiffmatWF_ADR1D_single(L,etas{i_dim},OO,OO,tau_grid(i_dim,:),opt);
    M_eta{i_dim}  = stiffmatWF_ADR1D_single(L,OO,OO,etas{i_dim},tau_grid(i_dim,:),opt);
    M_rho{i_dim}  = stiffmatWF_ADR1D_single(L,OO,OO,rhos{i_dim},tau_grid(i_dim,:),opt);
    T_beta{i_dim} = stiffmatWF_ADR1D_single(L,OO,betas{i_dim,i_dim},OO,tau_grid(i_dim,:),opt);
end

M_beta{1,1} = stiffmatWF_ADR1D_single(L,OO,OO,betas{1,2},tau_grid(2,:),opt);
M_beta{1,2} = stiffmatWF_ADR1D_single(L,OO,OO,betas{1,3},tau_grid(3,:),opt);
M_beta{2,1} = stiffmatWF_ADR1D_single(L,OO,OO,betas{2,1},tau_grid(1,:),opt);
M_beta{2,2} = stiffmatWF_ADR1D_single(L,OO,OO,betas{2,3},tau_grid(3,:),opt);
M_beta{3,1} = stiffmatWF_ADR1D_single(L,OO,OO,betas{3,1},tau_grid(1,:),opt);
M_beta{3,2} = stiffmatWF_ADR1D_single(L,OO,OO,betas{3,2},tau_grid(2,:),opt);



% assemble the 3D stiffness matrix
B = kron(kron(K_eta{1}, M_eta{2}), M_eta{3}) + ...
    kron(kron(M_eta{1}, K_eta{2}), M_eta{3}) + ...
    kron(kron(M_eta{1}, M_eta{2}), K_eta{3}) + ...
    kron(kron(T_beta{1},   M_beta{1,1}), M_beta{1,2}) + ...
    kron(kron(M_beta{2,1}, T_beta{2}),   M_beta{2,2}) + ...
    kron(kron(M_beta{3,1}, M_beta{3,2}), T_beta{3}) + ...
    kron(kron(M_rho{1}, M_rho{2}), M_rho{3});

end





