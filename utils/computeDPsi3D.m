function DPsi = computeDPsi3D(l0,L,type)
% DPSI = computeDPsi3D(L0,L,TYPE) computes the matrix associated with norm
% rescaling for wavelets of anisotropic or isotropic type. DPSI is a
% diagonal matrix whose entries are the inverses of the H1-norm of the
% wavelet functions

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

N1 = 2^L;
N  = N1^3;

switch type
    case 'ani'
        ll = max(l0,floor(log2(0:N1-1)));
        D_1D = sparse(diag(2.^ll));
        D_ani_norm = sparse(maxkron(maxkron(D_1D,D_1D),D_1D)); % |psi_{l,k}|_H1 ~ 2^norm(l,Inf)
        DPsi = sparse(diag(1./diag(D_ani_norm)));
    case 'iso'
        D = zeros(N1,N1,N1);
        for l = L:-1:l0
            D(1:2^l,1:2^l,1:2^l) = 2^(max(l-1,l0)); % |psi_{l,k}|_H1 ~ 2^l
        end
        DPsi = sparse(diag(1./D(:)));
end