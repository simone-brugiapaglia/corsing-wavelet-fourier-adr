%% Figure 7
% Figure 7A: Plot of the exact solution to the one-dimensional
% advection-diffusion-reaction problem and of its best s-term approximation
% Figure 7B: Plot of the wavelet coefficients of the exact solution and of
% its best s-term approximation

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

clear all
close all

addpath 'utils'

% exact solution
u_ex = @(x)  1+exp(-(x-0.3).^2/0.0005) + 0.5*(cos(2*pi*x));

% Trial space parameters
l0 = 2;
L = 9;
L_add = 3;
N = 2^L;
s = 50;
[primal,dual] = filterbank('bior'); % wavelet filters


%% Compute wavelet coefficients 
h = 2^(-L);
x = 0:h:1-h;
uu_scal = 2^(-L/2) * u_ex(x)'; % factor 2^(-L/2) is due to the fact that max(phi_Lk) = 2^(L/2)
uu_wave = wavedeco1D(uu_scal',l0,L,dual)';
ll = max(l0,floor(log2(0:N-1)));
DPsi = diag(2.^(-ll));
uu = uu_wave ./ diag(DPsi); % normalization wrt H1 norm


%% Compute best 50-term approximation error 
[uu_sorted,sorting] = sort(abs(uu),'descend');
uu_s = zeros(N,1);
uu_s(sorting(1:s)) = uu(sorting(1:s));
fprintf('Best %d-term approximation error = %1.4e\n',s,norm(uu - uu_s,2)/norm(uu,2))


%% Figure 1A
figure;
x = 0:h:1-h;
ux = u_ex(x);
usx = 2^(L/2)*wavereco1D((DPsi*uu_s)',l0,L,primal);
plot(x,ux,'--',x,usx,'-');
hl=legend('Exact solution',['Best ',num2str(s),'-term approximation']);
set(hl,'interpreter','latex')
pbaspect([2 1 1])
set(gca,'ticklabelinterpreter','latex')
set(gca,'fontsize',15)
grid on;
xlabel('$x$','interpreter','latex')
ylabel('$u(x)$','interpreter','latex')
ylim([0.4,2])
saveas(gca,['figures/',mfilename,'A'],'epsc')
    

%% Figure 1B
figure;
stem(uu); 
hold on; 
stem(sorting(1:s),uu_s(sorting(1:s)),'x');
hold off
hl=legend('Wavelet coefficients',['Best ',num2str(s),'-term approximation'],'location','northeast');
set(hl,'interpreter','latex')
grid on
set(gca,'ticklabelinterpreter','latex')  
set(gca,'fontsize',15)
axis tight
pbaspect([2 1 1])
xlabel('Index $j$','interpreter','latex')
ylabel('$j^{\mathrm{th}}$ entry','interpreter','latex')
saveas(gca,['figures/',mfilename,'B'],'epsc')
    