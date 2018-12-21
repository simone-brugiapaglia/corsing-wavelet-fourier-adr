%% Figure 9
% Sparsity vs. recovery error for a one-dimensional
% advection-diffusion-reaction probelm with constant and nonconstant
% coefficients.
%
% This figure is produced from precomputed data stored in 
% "/data/Data_Figure09_article". To repeat the computation, run the file
% "Compute_Figure09". This will produce a new MAT file
% "/data/Data_Figure09" that has to be loaded in place of
% "/data/Data_Figure09_article".

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)
clear all
close all

load data/Figure09_article.mat
%% Uncomment the following line to use new data
% load data/Figure08_new

%% Create figures
for i_eta = 1:2
    
    figure(i_eta);
    boxplot(squeeze(rel_err(i_eta,:,:))','labels',ss)
    set(gca,'ticklabelinterpreter','latex')
    grid on
    hold on
    plot(median(squeeze(rel_err(i_eta,:,:)),2),'-r');
    plot(best_s_term_err,'--');
    hold off
    hl = legend('CORSING','best $s$-term');
    set(hl,'interpreter','latex');
    xlabel('$s$','interpreter','latex');    
    ylim([min(best_s_term_err),max(max(max(rel_err)))])
    switch i_eta
        case 1
            ylabel('Relative $H^1(\mathcal{D})$-error','interpreter','latex')
    end
    set(gca, 'YScale', 'log')
    pbaspect([2 1 1])
    set(gca,'fontsize',20)
    switch i_eta
        case 1
            saveas(gca,['figures/',mfilename,'A'],'epsc')
        case 2
            saveas(gca,['figures/',mfilename,'B'],'epsc')
    end       
end