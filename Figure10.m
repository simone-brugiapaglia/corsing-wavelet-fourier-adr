%% Figure 10
% Visualization of the two-dimensional functions u2 and u3 with respective
% wavelet coefficients (anisotropic and isotropic).

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

clear all
close all

addpath utils

%% Low-rank representation of u1
u11 = { @(x)  exp(-(x-0.3).^2/0.0005), @(x) 2*exp(-(x-0.6).^2/0.001)}; 
u12 = { @(x)  exp(-(x-0.4).^2/0.0005), @(x) exp(-(x-0.5).^2/0.005)}; 
u{1} = @(x,y) lowrankfun2D(u11,u12,x,y);

%% Low-rank representation of u2
u21 = {@(x)  exp(-(x-0.45).^2/0.005)};
u22 = {@(x) 1+0*x};
u{2} = @(x,y) lowrankfun2D(u21,u22,x,y);

%% Evaluate exact solution at high-resolution uniform grid points 
h_hr= 0.001;
xh_hr = 0:h_hr:1;
[Xh1_hr,Xh2_hr] = meshgrid(xh_hr,xh_hr);

l0 = 2;
L = 6;
N1 = 2^L;
N = N1^2; % PG trial space dimension
n =2;
s = 100;
h = 2^(-L);
xh = 0:h:1-h;
[Xh1,Xh2] = meshgrid(xh,xh);

types = {'ani','iso'};
opt.trial = 'bior';
[primal,dual] = filterbank(opt.trial);


for i_u = 1:2 
    %% Function 2D plot
    figure;
    Uh_hr = u{i_u}(Xh1_hr,Xh2_hr);
    imagesc(xh_hr,xh_hr,Uh_hr);
    set(gca,'YDir','normal')
    set(gca,'ticklabelinterpreter','latex')
    hc = colorbar;
    set(hc,'ticklabelinterpreter','latex');
    set(gca,'fontsize',25)
    pbaspect([1 1 1])
    
    %% Save Figure
    switch i_u
        case 1
            saveas(gca,['figures/',mfilename,'A'],'epsc')
        case 2
            saveas(gca,['figures/',mfilename,'D'],'epsc')
    end
    
    %% Wavelet coefficients
    for i_type = 1:2
        Uh = u{i_u}(Xh1,Xh2);
        UU_single = (2^(-L/2))^n * Uh; % coefficients wrt the single-scale basis
        
        DPsi = computeDPsi2D(l0,L,types{i_type}); % DPsi ~ diag(1 / ||psi_j||_H1 )
        UU_wave = wavedeco2D(UU_single,l0,L,dual,types{i_type});     
        uu = reshape(UU_wave,N,1) ./ diag(DPsi); % normalization wrt H1 norm
        UU = reshape(uu,N1,N1);
        
        %% Compute best s-term approximation error
        uu_s = zeros(N,1);
        [u_sorted,sorting] = sort(abs(uu),'descend');
        uu_s(sorting(1:s)) = uu(sorting(1:s));
        
        rel_best_s_term = norm(uu - uu_s, 2) / norm(uu, 2);
        
        fprintf(['u',num2str(i_u),': best s-term (',types{i_type},') =  %1.4e\n'], rel_best_s_term);
        
        %% Plot wavelet coefficients
        figure;
        imagesc(abs(UU));
        set(gca,'YDir','normal')
        hc = colorbar;
        set(gca,'ticklabelinterpreter','latex');
        set(hc,'ticklabelinterpreter','latex');
        set(gca,'fontsize',25)
        pbaspect([1 1 1])
        
        %% Save Figure
        switch i_u
            case 1
                switch i_type
                    case 1
                        saveas(gca,['figures/',mfilename,'B'],'epsc')
                    case 2
                        saveas(gca,['figures/',mfilename,'C'],'epsc')
                end
            case 2
                switch i_type
                    case 1
                        saveas(gca,['figures/',mfilename,'E'],'epsc')
                    case 2
                        saveas(gca,['figures/',mfilename,'F'],'epsc')
                end
        end
                
    end
end
