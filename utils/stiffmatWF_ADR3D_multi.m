function B = stiffmatWF_ADR3D_multi(l0,L,etas,betas,rhos,tau_grid,opt)
%stiffmatWF_ADR3D_iso builds the Petrov-Galerkin stiffness matrix for the
% ADR bilinear form
%
% a(u,v) = (eta * grad(u), grad(v)) + (beta * grad(u) , v) + (rho * u , v),
%
% where the coefficients are separable 
%  eta(x,y)  = eta_1(x) eta_2(y) eta_3(z),
%  beta(x,y) = [beta_11(x), beta_12(y), beta_13(z) ; 
%               beta_21(x), beta_22(y), beta_23(z) ;
%               beta_31(x), beta_32(y), beta_33(z) ],
%  rho(x,y)  = rho_1(x) rho_2(y) rho_3(z).
% L     = max level
% etas  = {eta1,  eta2, eta3}
% betas = {beta1, beta2, beta3}
% rhos  = {rho1,  rho2}
% tau_grid   = the rows are the values of the grid of tests


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


%% Single-scale stiffness matrix
B_single = stiffmatWF_ADR3D_single(L,etas,betas,rhos,tau_grid,opt);

m = size(B_single,1);
N1 = 2^L;
N  = N1^3;

%% Single-scale to multi-scale conversion
[primal,dual] = filterbank(opt.trial);
B1 = zeros(size(B_single));
for i = 1:m
    CURR_ROW = reshape(B_single(i,:),N1,N1,N1);
    B1(i,:) = reshape(wavedeco3D(CURR_ROW,l0,L,primal,opt.type),1,N);
end

%% H1 normalization
[TAU1 , TAU2, TAU3] = meshgrid(tau_grid(1,:),tau_grid(2,:),tau_grid(3,:));
tau1 = TAU1(:);
tau2 = TAU2(:);
tau3 = TAU3(:);
DXi = sparse(diag(1./sqrt(1 + (2*pi)^2*(tau1.^2 + tau2.^2 + tau3.^2)))); % DXi ~ diag(1 / ||xi_j||_H1 )
DPsi = computeDPsi3D(l0,L,opt.type); % DPsi ~ diag(1 / ||psi_j||_H1 )
B = DXi * B1 * DPsi;


end

