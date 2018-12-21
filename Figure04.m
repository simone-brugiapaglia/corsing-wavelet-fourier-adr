%% Figure 4 
% This figure show the sharpness of the upper bound to the quantity
% |(D(phi_lk) , D(xi_q))| employed in the one-dimensional local a-coherence
% estimate.
%

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

clear all
close all

l = 5;
k = 4;

inner_prod = @(l,k,q) 4* 2.^(3*l/2) .* exp(-2i*pi*q.*k.*(2.^(-l))).*...
    sin(pi*q.*2.^(-l)).^2;

upper_bound = @(l,q,gam) 4* 2.^((3/2-gam).*l).* abs(pi*q).^gam;

q_min = -100;
q_max = 100;

qs = q_min : q_max;

%% Plot figure
plot(qs,max(abs(inner_prod(l,k,qs)),1e-1),'r')
hold on
for gam = 0:0.25:2
    plot(qs,upper_bound(l,qs,gam),'color',(gam+0.5)/2.5*[0 0 1]);
end
hold off

text(13,1600,'$\gamma = 2$',...
    'horizontalalignment','right',...
    'verticalalignment','top','color',(2+0.5)/2.5*[0 0 1],...
    'interpreter','latex','fontsize',15);

text(70,upper_bound(l,q_max,0),'$\gamma = 0$',...
    'horizontalalignment','right',...
    'verticalalignment','bottom','interpreter','latex','fontsize',15);

%% Formatting and style
hl = legend('$|(\varphi_{\ell,k}'',\xi_q'')|$','upper bound','location','northwest');
set(hl,'interpreter','latex')
xlabel('Index $q$','interpreter','latex')
ylim([0,2000])
pbaspect([2 1 1])
set(gca,'ticklabelinterpreter','latex')
set(gca,'fontsize',15)
grid on

%% Save figure
saveas(gca,['figures/',mfilename],'epsc')
