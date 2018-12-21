function C = maxkron(A,B)
%MAXKRON Summary of this function goes here
%   Detailed explanation goes here
%
% maxkron(A,B) is like a kron(A,B), but the product is replaced by the
% maximum operation
%
%      [                     ]
% C =  [    max{A_ij, B}     ]
%      [                     ]

if isequal(A,diag(diag(A))) && isequal(B,diag(diag(B)))
    m = size(A,1);
    d = [];
    for i = 1:m
        d = [d;max(A(i,i),diag(B))];
    end
    C = diag(d);
    
else
    error('Not implemented for nondiagonal matrices')
end

