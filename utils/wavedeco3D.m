function D = wavedeco3D(C,l0,L,filters,type)
% D = WAVEDECO3D(C,L0,L,FILTERS,TYPE) computes the wavelet decomposition of
% the three-dimensional array C with respect to the wavelet filters of
% anisotropic and isotropic type.
%
% INPUT
%  C        single-scale coefficients (3D-tensor)
%  FILTERS  wavelet filters as output by filterbank
%  TYPE     can be 'ani' or 'iso'.
%
% OUTPUT
%
%  D        multi-scale coefficients

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

if size(C,1) ~= size(C,2) || size(C,2) ~= size(C,3)
    error('Coefficient tensor should have equal sizes')
end

N1 = size(C,1);

switch type
    case 'ani'
        % Initialization
        D1   = zeros(N1,N1,N1);
        D12  = zeros(N1,N1,N1);
        D123 = zeros(N1,N1,N1);
        
        % Recall: wavedeco1D transforms a stack of rows
        
        % 1D transform across the first dimension
        for i3 = 1:N1
            D1(:,:,i3) = wavedeco1D(squeeze(C(:,:,i3))',l0,L,filters)';
        end
        
        % 1D transform across the second dimension
        for i3 = 1:N1
            D12(:,:,i3) = wavedeco1D(squeeze(D1(:,:,i3)),l0,L,filters);
        end
        
        % 1D transform across the third dimension
        for i1 = 1:N1
            D123(i1,:,:) = wavedeco1D(squeeze(D12(i1,:,:)),l0,L,filters);
        end
        
        D = D123;
    case 'iso'
        D = C;
        
        for l = L:-1:l0+1
            % extract single-scale decomposition al level l
            Cl = D(1:2^l,1:2^l,1:2^l);
            
            % initialization
            Dl1   = zeros(2^l,2^l,2^l);
            Dl12  = zeros(2^l,2^l,2^l);
            Dl123 = zeros(2^l,2^l,2^l);
            
            % 1D transform across the first dimension
            for i3 = 1:2^l
                Dl1(:,:,i3) = wavedeco1D(squeeze(Cl(:,:,i3))',l-1,l,filters)';
            end
            
            % 1D transform across the second dimension
            for i3 = 1:2^l
                Dl12(:,:,i3) = wavedeco1D(squeeze(Dl1(:,:,i3)),l-1,l,filters);
            end
            
            % 1D transform across the third dimension
            for i1 = 1:2^l
                Dl123(i1,:,:) = wavedeco1D(squeeze(Dl12(i1,:,:)),l-1,l,filters);
            end
            
            % substitute scaling coefficients at level l with scaling at level l-1 and details at level l
            D(1:2^l,1:2^l,1:2^l) = Dl123;
        end
        
        
end
        
