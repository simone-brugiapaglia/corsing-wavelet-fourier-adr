function gg = loadvecF3D(f1,f2,f3,tau_grid,h,NGpoints,norm)
% FF = LOADVECF3D(F,TAU,HINT,NGPOINTS) computes the load vector FF 
% with respect to the Fourier basis given the forcing term F. Numerical 
% integration is performed using a composite Gaussian formula of step H and
% with NGPOINTS Gaussian points.

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

%% Build normalization matrix
[TAU1 , TAU2, TAU3] = meshgrid(tau_grid(1,:),tau_grid(2,:),tau_grid(3,:));
tau1 = TAU1(:);
tau2 = TAU2(:);
tau3 = TAU3(:);


%% Assemble load vector
gg = zeros(size(tau_grid,2)^3,1);
for i_deco = 1:length(f1)
    gg_1{i_deco} = loadvecF1D(f1{i_deco},tau_grid(1,:),h,NGpoints,'L2');
    gg_2{i_deco} = loadvecF1D(f2{i_deco},tau_grid(2,:),h,NGpoints,'L2');
    gg_3{i_deco} = loadvecF1D(f3{i_deco},tau_grid(3,:),h,NGpoints,'L2');
    gg = gg + kron(kron(gg_1{i_deco},gg_2{i_deco}),gg_3{i_deco});
end
switch norm
    case 'H1'
        DXi = sparse(diag(1./sqrt(1 + (2*pi)^2*(tau1.^2 + tau2.^2 + tau3.^2)))); % DXi ~ diag(1 / ||xi_j||_H1 )
        gg  = DXi * gg;
    case 'L2'
    otherwise
        error('Norm type not valid.')
end