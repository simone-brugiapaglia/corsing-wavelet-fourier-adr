%% Figure 12
% Number of test vs. recovery error for CORSING WF applied to a
% two-dimensional adveciton-diffusion-reaction problem. Anisotropic and
% isotropic wavelet are considered, as well as uniform and nonuniform
% subsampling. The exact solution is anisotropic and is gaussian along the
% x1-direction and constant along the x2-direction.
%
% This figure is produced from precomputed data stored in 
% "/data/Data_Figure12_article". To repeat the computation, run the file
% "Compute_Figure12". This will produce a new MAT file
% "/data/Data_Figure12" that has to be loaded in place of
% "/data/Data_Figure12_article".

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)
clear all
close all

load data/Figure12_article.mat
%% Uncomment the following line to use new data
% load data/Figure11_new

i_fig = 0;
for i_trial = 1:length(trials)
    for i_type = 1:length(types)
        for i_nu = 1:length(nus)
            i_fig = i_fig + 1;
            figure(i_fig)
            boxplot(squeeze(rel_err(i_trial,i_type,i_nu,:,:))','labels',ms)
            set(gca,'ticklabelinterpreter','latex')
            grid on
            xlabel('$m$','interpreter','latex')
            if i_nu == 1 && i_type == 1
                ylabel('Relative $H^1(\mathcal{D})$-error','interpreter','latex')
            end
            ylim([min(min(min(min(min(rel_err))))),max(max(max(max(max(rel_err)))))])
            set(gca, 'YScale', 'log')
            pbaspect([1 1 1])
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
            if isequal([i_type,i_nu],[1,1])
                saveas(gca,['figures/',mfilename,'A'],'epsc')
            elseif isequal([i_type,i_nu],[1,2])
                saveas(gca,['figures/',mfilename,'B'],'epsc')
            elseif isequal([i_type,i_nu],[2,1])
                saveas(gca,['figures/',mfilename,'C'],'epsc')
            elseif isequal([i_type,i_nu],[2,2])
                saveas(gca,['figures/',mfilename,'D'],'epsc')
            end
            
        end
    end
end
