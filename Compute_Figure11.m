%% Computations for Figure 11
% Number of test vs. recovery error for CORSING WF applied to a
% two-dimensional adveciton-diffusion-reaction problem. Anisotropic and
% isotropic wavelet are considered, as well as uniform and nonuniform
% subsampling. The exact solution is a superposition of gaussians.
%
% Running this file will save the data in "data/Figure11_new". This new
% data can be plotted by modifying the script Figure11.m. 

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)
clear all;
close all;

addpath utils

figure_name = 'Figure11';

n = 2; % physical domain dimension

% Trial and test space parameters
trials = {'bior'};
types = {'ani','iso'};
l0 = 2;
L = 6; % 
N1 = 2^L;
N = N1^2; % PG trial space dimension
R = N1;
M = R^2; % PG test space dimension

%% CORSING parameters
s = 100;          % sparsity 
ms = 100:100:500;  % values for sampling complexity m
N_run = 100;      % number of random CORSING runs

% numerical integration parameters
h = 2^(-L);
opt.L_add = 0; % additional levels added to build A via wavedec
hint = h;          % integration stepsize (ff)
NGpoints = 2;      % Gauss point

%% PDE coefficients (diffusion, advection, reaction) assumed to be separable
etas = {@(x) 0*x + 1, @(y) 0*y + 1}; % eta = etas{1} x etas{2}
rhos = {@(x) 0*x + 1, @(y) 0*y + 1}; % rho = rhos{1} x rhos{2}
betas{1,1} = @(x) 0*x + 1;
betas{1,2} = @(y) 0*y + 1;
betas{2,1} = @(x) 0*x + 1;
betas{2,2} = @(y) 0*y + 1; % beta = [betas{1,1} x betas{1,2} ; betas{2,1} x betas{2,2}]

%% exact solution
u1 = { @(x)  exp(-(x-0.3).^2/0.0005), @(x) 2*exp(-(x-0.6).^2/0.001)}; 
u2 = { @(x)  exp(-(x-0.4).^2/0.0005), @(x) exp(-(x-0.5).^2/0.005)};  
u = @(x,y) lowrankfun2D(u1,u2,x,y);

%% CORSING probability distributions
QQ_grid_1 = floor(-R/2)+1:floor(R/2);
QQ_grid_2 = floor(-R/2)+1:floor(R/2);
[QQ1, QQ2] = meshgrid(QQ_grid_1,QQ_grid_2);
qq1 = QQ1(:);
qq2 = QQ2(:);
norm2_qq = sqrt(qq1.^2 + qq2.^2);
normInf_qq = max(abs(qq1), abs(qq2));
norm0_qq = (qq1 ~= 0) + (qq2 ~= 0);
prod_abs_qq = max(abs(qq1),1) .* max(abs(qq2),1);

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
[f1,f2] = forcingterm2D(u1,u2,etas,betas,rhos); 

% Evaluate exact solution at uniform grid points
xh = 0:h:1-h;
[Xh1,Xh2] = meshgrid(xh,xh);
Uh = u(Xh1,Xh2);

% for loop on trial families and types
for i_trial = 1:length(trials)
    opt.trial = trials{i_trial};
    [primal,dual] = filterbank(opt.trial);   % Wavelet filters
 
    
    for i_type = 1:length(types)
        disp([trials{i_trial},' ',types{i_type}]) 
        
        type = types{i_type}; 
        opt.type  = type; 
                       
        %% Precomputation: full Petrov-Galerkin WF stiffness matrix        
        B = stiffmatWF_ADR2D_multi(l0,L,etas,betas,rhos,[QQ_grid_1;QQ_grid_2],opt);       
        gg = loadvecF2D(f1,f2,[QQ_grid_1;QQ_grid_2],h,NGpoints,'H1');


        %% Compute wavelet coefficients
        DPsi = computeDPsi2D(l0,L,type);
        UU_single = (2^(-L/2))^n * Uh;
        UU_wave = wavedeco2D(UU_single,l0,L,dual,type);  
        uu = reshape(UU_wave,N,1) ./ diag(DPsi);

        
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
    
                    % Wavelet reconstruction                    
                    UU_CORSING_Xh = (2^(L/2))^n * wavereco2D(reshape(DPsi*uu_CORSING,N1,N1),l0,L,primal,type);
                    UU_CORSING = reshape(uu_CORSING,N1,N1);
                    
                    
                    % store relative error wrt the H1 norm
                    rel_err(i_trial,i_type,i_nu,i_m,i_run) = norm(uu_CORSING - uu, 2) / norm(uu, 2);
                end
            end
            
            %% storage
            UU_CORSING_Xhs{i_trial,i_type,i_nu} = UU_CORSING_Xh;
            uu_CORSINGs{i_trial,i_type,i_nu} = uu_CORSING;

               
            fprintf('\n')
        end
    end
end


save(['data/',figure_name,'_new'],'s','L','l0','N1','N','R','M','ms','rel_err','trials',...
    'types','nus','pdists','UU_CORSING_Xhs','uu_CORSINGs','Uh','Xh1','Xh2','opt',...
    'hint','NGpoints')



