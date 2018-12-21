function D = wavedeco2D(C,l0,L,filters,type)
% D = WAVEDECO2D(C,L0,L,FILTERS,TYPE) computes the wavelet decomposition of
% the two-dimensional array C with respect to the wavelet filters of
% anisotropic and isotropic type.
%
% INPUT
%  C        single-scale coefficients (2D-tensor)
%  FILTERS  wavelet filters as output by filterbank
%  TYPE     can be 'ani' or 'iso'.
%
% OUTPUT
%
%  D        multi-scale coefficients

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

switch type
    case 'ani'
        %% Anisotropic case
        % first wavelet decomposition by rows
        Drow = wavedeco1D(C,l0,L,filters);
        
        % further wavelet decomposition by columns
        D = wavedeco1D(Drow',l0,L,filters)';
        
    case 'iso'
        %% Isotropic case
        D = C;
        
        for l = L:-1:l0+1
            % extract multi-scale decomposition al level l
            C_l = D(1:2^l,1:2^l);
            
            % first wavelet decomposition by rows
            Drow_l = wavedeco1D(C_l,l-1,l,filters);
            
            % further wavelet decomposition by columns
            D_l = wavedeco1D(Drow_l',l-1,l,filters)';
            
            % substitute single scale with multi-scale at level l
            D(1:2^l,1:2^l) = D_l;
        end
    otherwise
        error('Wavelet type is not valid.')
end

