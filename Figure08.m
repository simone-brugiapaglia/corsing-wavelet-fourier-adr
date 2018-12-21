%% Figure 8
% Number of test functions vs. recovery error for a one-dimensional
% advection-diffusion-reaction probelm with constant and nonconstant
% coefficients.
%
% This figure is produced from precomputed data stored in 
% "/data/Data_Figure07_article". To repeat the computation, run the file
% "Compute_Figure08". This will produce a new MAT file
% "/data/Data_Figure08" that has to be loaded in place of
% "/data/Data_Figure08_article".

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

clear all
close all

load data/Figure08_article.mat
%% Uncomment the following line to use new data
% load data/Figure07_new

%% Figures
for i_eta = 1:2
    
    figure(i_eta);
    boxplot(squeeze(rel_err(i_eta,:,:))','labels',ms)
    hold on
    plot(0:length(ms)+1,1.68e-2*ones(1,length(ms)+2),'--') % comparison with best s-term approximation error
    hold off
    set(gca,'ticklabelinterpreter','latex')
    grid on
    xlabel('$m$','interpreter','latex')
    if i_eta == 1
        ylabel('Relative $H^1(\mathcal{D})$-error','interpreter','latex')
    end
    ylim([min([min(min(min(rel_err))),1.5e-2]),max(max(max(rel_err)))])

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