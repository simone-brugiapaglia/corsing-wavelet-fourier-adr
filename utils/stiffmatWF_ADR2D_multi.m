function B = stiffmatWF_ADR2D_multi(l0,L,etas,betas,rhos,tau_grid,opt)
%stiffmatWF_ADR2D_ani buils the Petrov-Galerkin stiffness matrix for the
% ADR bilinear form
%
% a(u,v) = (eta grad u, grad v) + (beta grad u , v) + (rho u , v),
%
% where the coefficients are separable 
%  eta(x,y)  = eta_1(x) eta_2(y),
%  beta(x,y) = [beta_11(x) beta_12(y) ; beta_21(x) beta_22(y)],
%  rho(x,y)  = rho_1(x) rho_2(y).
%
%
% l0         coarsest level
% L          max level
% etas{i}    eta_i
% betas{i,j} beta_ij
% rhos{i}    rho_i
% tau_grid  the rows are the values of the grid of tests
% opt   .forcequad
%       .L_add
%       .trial
%       .type

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

if nargin == 6
    opt.forcequad = 0;
    opt.L_add = 0;
    opt.trial = 'bior';
    opt.type = 'ani';
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
    if ~isfield(opt,'type')
        opt.type = 'ani';
    end
end

if (~isconst(etas{1}) || ~isconst(etas{2}) || ~isconst(betas{1,1}) || ...
        ~isconst(betas{1,2}) || ~isconst(betas{2,1}) || ...
        ~isconst(betas{2,2}) || ~isconst(rhos{1}) || ~isconst(rhos{2}) )
    error('2D implementation tensted only for constant coefficients')
end



%% Single-scale stiffness matrix
B_single = stiffmatWF_ADR2D_single(L,etas,betas,rhos,tau_grid,opt);

m = size(B_single,1);
N1 = 2^L;
N2 = N1^2;

%% Single-scale to multi-scale conversion
[primal,dual] = filterbank(opt.trial);
B1 = zeros(size(B_single));
for i = 1:m
    CURR_ROW = reshape(B_single(i,:),N1,N1);
    B1(i,:) = reshape(wavedeco2D(CURR_ROW,l0,L,primal,opt.type),1,N2);
end

%% Normalize w.r.t. the H1 norm
DPsi = computeDPsi2D(l0,L,opt.type); % DPsi ~ diag(1 / ||psi_j||_H1 )

[TAU1, TAU2] = meshgrid(tau_grid(1,:),tau_grid(2,:));
tau1 = TAU1(:);
tau2 = TAU2(:);
DXi  = sparse(diag(1./sqrt(1 + (2*pi)^2*(tau1.^2+tau2.^2)))); % DXi ~ diag(1 / ||xi_j||_H1 )

B = DXi * B1 * DPsi;
