function ff = loadvecF1D(f,tau,h,NGpoints,norm)
% FF = LOADVECF1D(F,TAU,HINT,NGPOINTS) computes the load vector FF 
% with respect to the Fourier basis given the forcing term F. Numerical 
% integration is performed using a composite Gaussian formula of step H and
% with NGPOINTS Gaussian points. The test functions are normalized with
% respect to the norm in NORM.
%
% F        function handle
% TAU      vector
% H        real scalar
% NGpoints integer
% NORM     string. Can be 'H1' or 'L2'.

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

if nargin < 5
    error('Not enough input arguments')
end

m = length(tau);

ff = zeros(m,1);
for i = 1:m
    phi  = @(x) exp(2i*pi*tau(i)*x);
    ff(i) = GaussInt1D(@(x) f(x).*conj(phi(x)),h,NGpoints);
    switch norm
        case 'H1'
            norm_phi = sqrt(1+(2*pi*tau(i))^2);
            ff(i) = ff(i) / norm_phi;
        case 'L2'
        otherwise
            error('Norm type not valid.')
    end
end
end