function C = wavereco3D(D,l0,L,filters,type)
% D = WAVERECO3D(C,L0,L,FILTERS,TYPE) computes the wavelet reconstruction
% of the three-dimensional array D with respect to the wavelet filters of
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


if size(D,1) ~= size(D,2) || size(D,2) ~= size(D,3)
    error('Coefficient tensor should have equal sizes')
end

N1 = size(D,1);

switch type
    case 'ani'
        
        C1   = zeros(N1,N1,N1);
        C12  = zeros(N1,N1,N1);
        C123 = zeros(N1,N1,N1);
        
        % Recall: wavereco1D transforms a stack of rows
        
        % 1D transform across the first dimension
        for i3 = 1:N1
            C1(:,:,i3) = wavereco1D(squeeze(D(:,:,i3))',l0,L,filters)';
        end
        
        % 1D transform across the second dimension
        for i3 = 1:N1
            C12(:,:,i3) = wavereco1D(squeeze(C1(:,:,i3)),l0,L,filters);
        end
        
        % 1D transform across the third dimension
        for i1 = 1:N1
            C123(i1,:,:) = wavereco1D(squeeze(C12(i1,:,:)),l0,L,filters);
        end
        
        C = C123;
        
    case 'iso'
        C = D;
        
        for l = l0:L-1
            % extract two-scale decomposition (scaling level l + details level l+1)
            Dl = C(1:2^(l+1),1:2^(l+1),1:2^(l+1));
            
            % initialization
            Cl1   = zeros(2^(l+1),2^(l+1),2^(l+1));
            Cl12  = zeros(2^(l+1),2^(l+1),2^(l+1));
            Cl123 = zeros(2^(l+1),2^(l+1),2^(l+1));
            
            % 1D transform across the first dimension
            for i3 = 1:2^(l+1)
                Cl1(:,:,i3) = wavereco1D(squeeze(Dl(:,:,i3))',l,l+1,filters)';
            end
            
            % 1D transform across the second dimension
            for i3 = 1:2^(l+1)
                Cl12(:,:,i3) = wavereco1D(squeeze(Cl1(:,:,i3)),l,l+1,filters);
            end
            
            % 1D transform across the second dimension
            for i1 = 1:2^(l+1)
                Cl123(i1,:,:) = wavereco1D(squeeze(Cl12(i1,:,:)),l,l+1,filters);
            end
            
            % substitute two-scale decomposition (scaling l-1 and details l) with single-scale at level l
            C(1:2^(l+1),1:2^(l+1),1:2^(l+1)) = Cl123;
        end
end
