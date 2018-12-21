function C = wavereco1D(X,l0,L,filters)
% WAVERECO computes the single scale coefficients from the vector
% associated with the multiscale decomposition x.
%
% INPUT:
%
%  x        multiscale coefficients
%  l0       coarsest level of the MRA
%  L        finest level of the MRA
%  filters  wavelet filters
%
% OUTPUT:
% 
%  D   each row contains the output of the cascade algorithm associated
%  with the corresponding row of X
%
%


% Simone Brugiapaglia, Fabio Nobile
% EPFL, 2016 - SFU, 2018

% Extract wavelet filters
a = filters.a;
b = filters.b;
a_start = filters.a_start;
b_start = filters.b_start;

la = length(a);
lb = length(b);

nrows = size(X,1);

C = X(:,1:2^l0); % single scale coefficients at the coarsest level l0  

% customized commands
rowfft    = @(X) fft(X,[],2);  % performs the fft of each row
rowifft   = @(X) ifft(X,[],2); % performs the ifft of each row
replicate = @(v) repmat(v,[nrows,1]); % replicates v on 'nrows' rows  
convo     = @(a,X) rowifft(conj(replicate(fft(a))).*rowfft(X));

for l = l0:L-1
    % extract wavelet coefficients at level l
    D = X(:,2^l+1:2^(l+1));
    
    % zero padding and shifting a 
    a_zp = zeros(1,2^(l+1));
    a_zp(1:la) = fliplr(a);
    a_zp = circshift(a_zp, [0 -la - a_start + 1 ]);
    
    % zero padding and shifting b
    b_zp = zeros(1,2^(l+1));
    b_zp(1:lb) = fliplr(b);
    b_zp = circshift(b_zp, [0 -lb - b_start + 1]);
    
    % upsampling
    C_up = zeros(nrows,2^(l+1));
    C_up(:,1:2:end) = C;
    
    D_up = zeros(nrows,2^(l+1));
    D_up(:,1:2:end) = D;
    
    % single transform step
    C = 1/sqrt(2) * (convo(a_zp,C_up) + convo(b_zp,D_up));

end

