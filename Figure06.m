%% Figure 6
% Visualization of the stiffness matrix associated with the wavelet-Fourier
% Petrov-Galerkin discretization of a one-dimensional
% advection-diffusion-reaction problem with constant and nonconstant
% diffusion coefficient.

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

clear all
close all

addpath 'utils'

%% PDE coefficients (diffusion, advection, reaction)
etas = {@(x) 1+0*x, @(x)  1 + 0.5*sin(6*pi*x)}; % constant and nonconstant diffusion coeficient
beta = @(x) 0*x ;
rho  = @(x) 0*x + 1 ;

% Trial space parameters
opt.trial = 'bior';
l0 = 2;
L = 9;
L_add = 3;
N = 2^L;
opt.L_add = L_add; % additional levels to increase integration accuracy

% Test space parameters
R = N;
qq = floor(-R/2)+1:floor(R/2); % test indices


for i_eta = 1:2
    eta = etas{i_eta};
    
    %% Stiffness matrix assembly
    tic
    B = stiffmatWF_ADR1D_multi(l0,L,eta,beta,rho,qq,opt); 
    t_stiffness = toc;
    
    DXi  = diag(sqrt(1 + (2i*pi*qq).^2));
    ll = max(l0,floor(log2(0:N-1)));
    DPsi = diag(2.^(ll));
    B_not_norm = DXi * B * DPsi;
    
    %% Print condition number info
    switch i_eta
        case 1
            fprintf('Constant diffusion:\n')
        case 2
            fprintf('Nonconstant diffusion:\n')
    end
    fprintf('cond(B) (without H1-normalization) = %1.3e\n',cond(B_not_norm))
    fprintf('cond(B) (with H1-normalization)    = %1.3e\n',cond(B))
    
    %% Visualize stiffness matrix
    figure
    imagesc(1:N, qq, abs(B)); 
    pbaspect([1 1 1])
    hc = colorbar;
    colormap parula
    switch i_eta
        case 1
            ht = title('$|B_{qj}|$');
        case 2
            ht = title('$|B_{qj}|$');
    end
    xlabel('$j$','interpreter','latex')
    ylabel('$q$','interpreter','latex')
    set(hc,'TickLabelInterpreter', 'latex');
    set(ht,'interpreter','latex');
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'YDir','normal')
    set(gca,'fontsize',15)
    
    %% Save figure
    switch i_eta
        case 1
            saveas(gca,['figures/',mfilename,'A'],'epsc')
        case 2
            saveas(gca,['figures/',mfilename,'B'],'epsc')
    end
    
end