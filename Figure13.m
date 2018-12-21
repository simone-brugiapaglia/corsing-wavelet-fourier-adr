%% Figure 13
% Visualization of the (anisotropic and isotropic) wavelet coefficients of 
% the three-dimensional functions used to test CORSING.

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

clear all
close all

addpath utils

l0 = 2;
L = 4;
s = 200;
n = 3; % dimension of the physical domain

N1 = 2^L;
N = N1^n;

%% Exact solutions
u1 = {@(x) exp(-(x-0.4).^2/0.005)};
u2 = {@(x) exp(-(x-0.5).^2/0.0005)};
u3 = {@(x) exp(-(x-0.6).^2/0.005)};
u{1} = @(x,y,z) lowrankfun3D(u1,u2,u3,x,y,z);


[primal,dual] = filterbank('bior');

%% Compute wavelet coefficients
types = {'ani','iso'};
h = 2^(-L);


for i_u = 1:length(u)
    xh = 0:h:1-h;
    [Xh,Yh,Zh] = meshgrid(xh,xh,xh);
    Uh = u{i_u}(Xh,Yh,Zh);
    UU_single = (2^(-L/2))^n * Uh;
    UU_single = permute(UU_single,[3,1,2]); % Adjust ordering to have accordance between
                                            % meshgrid and functions implemented in utils 
    for i_type = 1: length(types)
        %% Compute wavelet coefficients
        type = types{i_type};       
        Uwave = wavedeco3D(UU_single,l0,L,dual,type);
        DPsi = computeDPsi3D(l0,L,type);
        uu = Uwave(:)./diag(DPsi); % H1-normalization
        
        %% compute best s-term
        uu_s = zeros(size(uu));
        [uu_sorted,sorting] = sort(abs(uu),'descend');
        uu_s(sorting(1:s)) = uu(sorting(1:s));
        rel_best_s_err = norm(uu-uu_s,2) / norm(uu,2);
        fprintf(['u',num2str(i_u),': rel best s-term (',type,') = %1.4e\n'],rel_best_s_err)
        
        %% Reconstruction
        UUs_reco = wavereco3D(reshape(DPsi*uu_s,N1,N1,N1),l0,L,primal,type);
        UUs_h = (2^(L/2))^n * permute(UUs_reco,[2,3,1]);
        
        
        fprintf('Norm infinity rel error = %1.4e\n',norm(Uh(:) - UUs_h(:), Inf) / norm(Uh(:),Inf));
        
        %% Plot coefficients
        figure;
        stem(uu);
        hold on
        stem(sorting(1:s),uu_s(sorting(1:s)),'x')
        hold off
        
        
      
        %% Visual settings
        grid on
        hl = legend('Wavelet coefficients',['Best ',num2str(s),'-term approx.'],'location','northeast');
        set(hl,'interpreter','latex')
        xlabel('Index $j$','interpreter','latex')

        if i_type == 1
            ylabel('$j^{\mathrm{th}}$ entry','interpreter','latex')
        end
        set(gca,'ticklabelinterpreter','latex')
        axis tight
        ylim([min(uu),max(uu)])
        pbaspect([2 1 1])
        set(gca,'fontsize',15)
        
        if isequal([i_u,i_type],[1,1])
            saveas(gca,['figures/',mfilename,'A'],'epsc')
        elseif isequal([i_u,i_type],[1,2])
            saveas(gca,['figures/',mfilename,'B'],'epsc')
        elseif isequal([i_u,i_type],[2,1])
            saveas(gca,['figures/',mfilename,'C'],'epsc')
        elseif isequal([i_u,i_type],[2,2])
            saveas(gca,['figures/',mfilename,'D'],'epsc')
        end
        
        %figure;
        %subplot(1,2,1);
        %slice(Xh,Yh,Zh,Uh,0.4,0.5,0.6);
        %colorbar;
        
        %subplot(1,2,2);
        %slice(Xh,Yh,Zh,UUs_h,0.5,0.6,0.4);
        %colorbar;
        
    end
end