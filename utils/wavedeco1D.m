function D = wavedeco1D(X,l0,L,filters)
% WAVEDECO converts the rows of X from single scale representation to 
% multiscale decomposition.
%
% INPUT:
%
%  X        each row contains the coefficients at the finest level L
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
% EPFL, 2016
   

% Extract wavelet filters
a = filters.a;
b = filters.b;
a_start = filters.a_start;
b_start = filters.b_start;

nrows = size(X,1);

la = length(a);
lb = length(b);

D = [];

% customized operations
rowfft    = @(X) fft(X,[],2);  % performs the fft of each row
rowifft   = @(X) ifft(X,[],2); % performs the ifft of each row
replicate = @(v) repmat(v,[nrows,1]); % replicates v on 'nrows' rows  

for l = L:-1:l0+1
    % zero padding a
    a_zp= zeros(1,2^l);
    a_zp(1:la) = a;
    a_zp = circshift(a_zp,[0 a_start]);
    
    % zero padding b
    b_zp = zeros(1,2^l);
    b_zp(1:lb) = b;
    b_zp = circshift(b_zp,[0 b_start]);
    
    Y = 1/sqrt(2) * rowifft(conj(replicate(fft(b_zp))).*rowfft(X)); % convolution
    D = [Y(:,1:2:end-1), D]; % downsampling y (notice that 0 in Matlab is 1)
    
    C = 1/sqrt(2) * rowifft(conj(replicate(fft(a_zp))).*rowfft(X)); % convolution
    C = C(:,1:2:end-1); % downsampling c (notice that 0 in Matlab is 1)
    
    X = C; % recursive step
end
D = [C, D]; % wavelet decomposition [c^l0, d^l0, ..., d^(L-1)] of ctilde

