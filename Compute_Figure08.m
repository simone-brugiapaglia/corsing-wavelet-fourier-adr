%% Computation of Figures 8A and 8B
% Comparison between CORSING on a problem with constant diffusion 
% coefficient and CORSING on a problem with nonconstant diffusion
% coefficient. The two cases are compared by studying the recovery error
% as a function of the numeber of tests m.
%
% Running this file will save the data in "data/Figure08_new". This new
% data can be plotted by modifying the script Figure08.m.

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

clear all
close all

addpath 'utils'
figure_name = 'Figure08';

%% Parameter setting
% ADR coefficients
etas = {@(x) 1+0*x, @(x)  1 + 0.5*sin(6*pi*x)};
beta = @(x) 0*x ;
rho  = @(x) 0*x + 1 ;

% exact solution
u = @(x)  1+exp(-(x-0.3).^2/0.0005) + 0.5*(cos(2*pi*x));

% trial functions settings
l0 = 2;
L = 9; 
L_add = 3;
N = 2^L;
opt.trial = 'bior';
[primal,dual] = filterbank('bior');
ll = max(l0,floor(log2(0:N-1)));
DPsi = diag(2.^(-ll)); % Diagonal matrix ~ 1/||psi_j||_H1

% Test space
R = N; % PG test space dimension (square PG matrix)
qq = floor(-R/2)+1:floor(R/2);

s = 50;
ms = 2*s:50:N;
N_runs = 100; % number of randomized CORSING runs 


% numerical integration parameters (for computing load vector)
h = 2^(-L);
hint = h;             % integration stepsize (ff)
NGpoints = 2;         % Gauss points
opt.L_add = L_add;    % additional levels added to build A via wavedec


%% compute wavelet coefficients by wavelet decomposition
x = 0:h:1-h;
uu_scal = 2^(-L/2) * u(x)'; % single scale coefficients as row (rescaling 2^(-L/2) because phi_{L,k}(kh)=2^(L/2)
uu_wave = wavedeco1D(uu_scal',l0,L,dual)';
uu = uu_wave ./ diag(DPsi); % normalization wrt H1 norm

% Compute best 50-term approximation error wrt H1 norm
[uu_sorted,sorting] = sort(abs(uu),'descend');
uu_s = zeros(N,1);
uu_s(sorting(1:s)) = uu(sorting(1:s));
best_s_term_err = norm(uu - uu_s,2)/norm(uu,2);

% CORSING probability distribution
nu = min(1./abs(qq),1);
pdist = nu/sum(nu);
cdist = cumsum(pdist);

rel_err = zeros(length(etas),length(ms),N_runs);

for i_eta = 1:2
    eta = etas{i_eta};   
    f = forcingterm1D(u,eta,beta,rho); % compute forcing term symbolically
    
    %% Precomputation of the full Petrov-Galerkin discretization 
    B  = stiffmatWF_ADR1D_multi(l0,L,eta,beta,rho,qq,opt); % stiffness matrix
    gg = loadvecF1D(f,qq,hint,NGpoints,'H1'); % load vector     
    
    switch i_eta
        case 1
            fprintf('Constant diffusion:\n m = ')
        case 2
            fprintf('Nonconstant diffusion:\n m = ')
    end
    
    %% For loop with multiple CORSING runs
    i_m =0;
    for m = ms
        i_m = i_m +1;
        fprintf('%d ',m);        
        for i_test = 1:N_runs            
            %% CORSING
            tau = nurandi(R,m,pdist,cdist); % randomized test selection
            D = diag(1./sqrt(m*pdist(tau))); % Diagonal precodnitioner
            A = B(tau,:); % stiffness matrix
            ff = gg(tau); % load vector          
            
            % recovery via OMP (N.B. ompbox needs real coefficients)
            A_aug  = [real(D*A); imag(D*A)]; % Augmented real matrix
            ff_aug = [real(D*ff); imag(D*ff)]; % Augmented real load vector
            [Anorm,norms] = normalize(A_aug,'col');
            uu_OMP = omp(Anorm'*ff_aug,Anorm'*Anorm,s);
            uu_CORSING = uu_OMP./norms;
            
            % store relative recovery error 
            rel_err(i_eta,i_m,i_test) = norm(uu_CORSING - uu ,2) / norm(uu ,2); 
        end
    end
    fprintf('\n')
end

save(['data/',figure_name,'_new'],'s','N','R','ms','rel_err','best_s_term_err')
