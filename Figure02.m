%% Figure 2 
% This figure shows the difference between anisotropic and isotropic tensor
% product wavelets
%

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

clear all
close all

addpath utils

FONTSIZE = 15;

l0 = 2;
L = 5;
D = zeros(2^L,2^L);
D(2,30) = 1;
[primal,dual] = filterbank('bior');

for type = [1,2]
    figure
    h = 2^(-L);
    xgrid = 0:h:1-h;
    switch type
        case 1
            C = wavereco2D(D,l0,L,primal,'ani');
            surf(xgrid,xgrid,C)
            ht = title('$\psi_{(4,2),(13,1)}^{\mathrm{ani}} = \psi_{4,13} \otimes \varphi_{2,1}$');
            set(ht,'interpreter','latex')
            
        case 2
            C = wavereco2D(D,l0,L,primal,'iso');
            surf(xgrid,xgrid,C)
            ht = title('$\psi_{4,(13,1),(1,0)}^{\mathrm{iso}} = \psi_{4,13} \otimes \varphi_{4,1}$');
            set(ht,'interpreter','latex')
    end
    
    
    view(-17,17)
    axis tight
    axis equal
    %colorbar
    colormap winter
    xlabel('$x_1$','interpreter','latex')
    ylabel('$x_2$','interpreter','latex')
    set(gca,'ticklabelinterpreter','latex')
    
    
    
    switch type
        case 1
            set(gca,'fontsize',FONTSIZE+3);
            saveas(gca,['figures/',mfilename,'A'],'epsc')
        case 2
            set(gca,'fontsize',FONTSIZE);
            saveas(gca,['figures/',mfilename,'B'],'epsc')
    end
   
end