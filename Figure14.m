%% Figure 14
% Number of test vs. recovery error for CORSING WF applied to a
% two-dimensional adveciton-diffusion-reaction problem. Anisotropic and
% isotropic wavelet are considered, as well as uniform and nonuniform
% subsampling. The exact solution is a superposition of gaussians.
%
% This figure is produced from precomputed data stored in 
% "/data/Data_Figure14_article". To repeat the computation, run the file
% "Compute_Figure14". This will produce a new MAT file
% "/data/Data_Figure14" that has to be loaded in place of
% "/data/Data_Figure14_article".

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)
clear all
close all

load data/Figure14_article.mat
%% Uncomment the following line to use new data
%load data/Figure13_new


for i_type = 1:length(types)
    for i_nu = 1:length(nus)
        figure;
        boxplot(squeeze(rel_err(i_type,i_nu,:,:))','labels',ms)
        ylim([min(min(min(min(rel_err)))),max(max(max(max(rel_err))))])
        
        set(gca, 'YScale', 'log')
        set(gca,'ticklabelinterpreter','latex')
        grid on
        xlabel('$m$','interpreter','latex')
        %ylabel('')
        pbaspect([1 1 1])
        if i_nu == 1 && i_type == 1
            ylabel('Relative $H^1(\mathcal{D})$-error','interpreter','latex')
           
        end
        if isequal([i_type,i_nu],[1,1])
            ht = title('Anisotropic -- Uniform $\mathbf{p}$');
        elseif isequal([i_type,i_nu],[1,2])
            ht = title('Anisotropic -- Nonuniform $\mathbf{p}$');
        elseif isequal([i_type,i_nu],[2,1])
            ht = title('Isotropic -- Uniform $\mathbf{p}$');
        elseif isequal([i_type,i_nu],[2,2])
            ht = title('Isotropic -- Nonuniform $\mathbf{p}$');
        end
        set(ht,'interpreter','latex')
        set(gca,'fontsize',25)
        
        SIZE = [18,18];
        set(gcf,'units','centimeters')
        set(gcf, 'PaperUnits', 'centimeters');        
        set(gcf, 'Position', [1,1,SIZE(1),SIZE(2)])        
        set(gcf, 'PaperSize', [SIZE(1) SIZE(2)]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 SIZE(2)/100 SIZE(1) 99/100*(SIZE(2))]);
        set(gcf, 'renderer', 'painters');
        
        
        if isequal([i_type,i_nu],[1,1])
            saveas(gcf,['figures/',mfilename,'A'],'epsc')
        elseif isequal([i_type,i_nu],[1,2])
            saveas(gcf,['figures/',mfilename,'B'],'epsc')
        elseif isequal([i_type,i_nu],[2,1])
            saveas(gcf,['figures/',mfilename,'C'],'epsc')
        elseif isequal([i_type,i_nu],[2,2])
            saveas(gcf,['figures/',mfilename,'D'],'epsc')
        end
        
    end
end

