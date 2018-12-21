function C = wavereco2D(D,l0,L,filters,type)
% D = WAVERECO2D(C,L0,L,FILTERS,TYPE) computes the wavelet reconstruction 
% of the two-dimensional array D with respect to the wavelet filters of
% anisotropic and isotropic type.

% INPUT
%  D        multi-scale coefficients
%  FILTERS  wavelet filters as output by filterbank
%  TYPE     can be 'ani' or 'iso'.
%
% OUTPUT
%  
%  C        single-scale coefficients

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)


switch type
    case 'ani'
        %% Anisotropic case
        % first wavelet decomposition by rows
        Crow = wavereco1D(D,l0,L,filters);
        
        % further wavelet decomposition by columns
        C = wavereco1D(Crow',l0,L,filters)';
        
    case 'iso' 
        %% Isotropic case
        C = D;
        
        for l = l0:L-1
            % extract multi-scale decomposition al level l
            D_l = C(1:2^(l+1),1:2^(l+1));
            
            % first wavelet decomposition by rows
            Crow_l = wavereco1D(D_l,l,l+1,filters);
            
            % further wavelet decomposition by columns
            C_l = wavereco1D(Crow_l',l,l+1,filters)';
            
            % substitute multi-scale with single scale at level l
            C(1:2^(l+1),1:2^(l+1)) = C_l;
        end
    otherwise
        error('Wavelet type is not valid.')
end