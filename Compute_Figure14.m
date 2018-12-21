%% Computations for Figure 14
% Number of test vs. recovery error for CORSING WF applied to a
% three-dimensional adveciton-diffusion-reaction problem. Anisotropic and
% isotropic wavelet are considered, as well as uniform and nonuniform
% subsampling. The exact solution is trigonometric.
%
% Running this file will save the data in "data/Figure14_new". This new
% data can be plotted by modifying the script Figure14.m. 

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)
clear all;
close all;

%% Trial space parameters
n = 3; % physical domain dimension
l0 = 2;
L = 4;
s = 200;
N1 = 2^L;
N = N1^3;
h = 2^(-L);
types = {'ani','iso'};
opt.trial = 'bior';
[primal,dual] = filterbank('bior');

N_run = 100;
ms = 200:100:600;

%% Test space parameters
R = N1;
M = R^3;

%% Exact solution
u1 = {@(x) cos(2*pi*x)};
u2 = {@(x) cos(2*pi*x)+sin(2*pi*x)};
u3 = {@(x) sin(2*pi*x)};
u = @(x,y,z) lowrankfun3D(u1,u2,u3,x,y,z);

%% Numerical integration parameters
hint = h;          % integration stepsize (ff)
NGpoints = 2;      % Gauss point


%% PDE coefficients (diffusion, advection, reaction) assumed to be separable
etas = {@(x) 0*x + 1, @(y) 0*y + 1, @(z) 0*z + 1}; % eta = etas{1} x etas{2}
rhos = {@(x) 0*x + 1, @(y) 0*y + 1, @(z) 0*z + 1}; % rho = rhos{1} x rhos{2}
for i = 1:3
    for j = 1:3
        betas{i,j} = @(x) 0*x + 1;
    end
end




%% CORSING probability distributions
QQ_grid_1 = floor(-R/2)+1 : floor(R/2);
QQ_grid_2 = floor(-R/2)+1 : floor(R/2);
QQ_grid_3 = floor(-R/2)+1 : floor(R/2);
[QQ1, QQ2, QQ3] = meshgrid(QQ_grid_1,QQ_grid_2,QQ_grid_3);
qq1 = QQ1(:);
qq2 = QQ2(:);
qq3 = QQ3(:);
norm2_qq = sqrt(qq1.^2 + qq2.^2 + qq3.^2);
normInf_qq = max(max(abs(qq1), abs(qq2)), abs(qq3));
norm0_qq = (qq1 ~= 0) + (qq2 ~= 0) + (qq3 ~= 0);
prod_abs_qq = max(abs(qq1),1) .* max(abs(qq2),1) .* max(abs(qq3),1);

DXi  = sparse(diag(1./sqrt(1 + (2*pi)^2 * (qq1.^2 + qq2.^2 + qq3.^2)))); % DXi ~ diag(||xi_q||_H1) 

% upper bound 1 (trial upper bound leading to the uniform probability)
nus{1} = ones(size(qq1));

% upper bound 2 (from the article, leading to nonuniform probability)
nus{2} = min(( norm2_qq.^2 ) ./ ( normInf_qq.^2 .* prod_abs_qq )  , 1 ); % local a-coherence upper bound

% compute probability and cumulative distributions
for i_nu = 1:length(nus)
    pdists{i_nu} = nus{i_nu}/sum(nus{i_nu}); % distribution
    cdists{i_nu} = cumsum(pdists{i_nu});
end


% Symbolic computation of the forcing term
[f1,f2,f3] = forcingterm3D(u1,u2,u3,etas,betas,rhos); 

for i_type = 1:length(types)
    type = types{i_type};
    opt.type = type;

    
    %% Compute wavelet coefficients 
    DPsi = computeDPsi3D(l0,L,type); % DPsi ~ diag(||psi_j||_H1)
    xh = 0:h:1-h;
    [Xh1,Xh2,Xh3] = meshgrid(xh,xh,xh);
    UU_single = (2^(-L/2))^n * u(Xh1,Xh2,Xh3); % Coefficients w.r.t. the single-scale basis
                                               % Rescaling is due because single-scale basis 
                                               % has max value 2^(L/2)^3
    UU_single = permute(UU_single,[3,1,2]);  % Adjust ordering to have accordance between
                                             % meshgrid and functions implemented in utils   
    UU_wave = wavedeco3D(UU_single,l0,L,dual,type);
    uu = reshape(UU_wave,N,1) ./ diag(DPsi); % normalization wrt H1 norm
    UU = reshape(uu,N1,N1,N1);   
    
    %% Precomputation: Stiffness matrix and load vector assembly
    B = stiffmatWF_ADR3D_multi(l0,L,etas,betas,rhos,[QQ_grid_1; QQ_grid_2; QQ_grid_3],opt);
    gg = loadvecF3D(f1,f2,f3,[QQ_grid_1; QQ_grid_2; QQ_grid_3],h,NGpoints,'H1');
    
        
    for i_nu = 1:length(nus)
        
        disp([ ' nu ',num2str(i_nu)]);
        
        pdist = pdists{i_nu};
        cdist = cdists{i_nu};
        
        
        %% For loop with multiple CORSING runs
        i_m =0;
        for m = ms
            i_m = i_m +1;
            fprintf('%d ',m)
            for i_run = 1:N_run
                %% CORSING
                tau = nurandi(N,m,pdist,cdist); % randomized test selection
                D = diag(1./sqrt(m*pdist(tau))); % Diagonal precodnitioner
                A = B(tau,:); % stiffness matrix
                ff = gg(tau); % load vector
                
                % recovery via OMP (N.B. ompbox needs real coefficients)
                A_aug  = [real(D*A);  imag(D*A) ]; % Augmented real matrix
                ff_aug = [real(D*ff); imag(D*ff)]; % Augmented real load vector
                [Anorm, norms] = normalize(A_aug,'col');
                uu_OMP = omp(Anorm'*ff_aug,Anorm'*Anorm,s);
                uu_CORSING = uu_OMP./norms;
                               
                % store relative error wrt the H1 norm
                rel_err(i_type,i_nu,i_m,i_run) = norm(uu_CORSING - uu, 2) / norm(uu, 2);
            end
        end        
        
        
        fprintf('\n')
    end
end


save('data/Figure14_new','s','L','l0','N1','N','R','M','ms','rel_err',...
    'types','nus','pdists','opt','hint','NGpoints','n')
