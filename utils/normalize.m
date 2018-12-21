function [B,norms] = normalize(A,opt)
% B = NORMALIZE(A,'row') normalizes the rows of A
%
% B = NORMALIZE(A,'col') normalizes the columns of A
%
% [B,NORMS] = NORMALIZE(A,OPT) returns the normalized matrix B and a
% vector NORMS storing the norms of the rows (or columns) of A.

% Simone Brugiapaglia, SFU
% simone_brugiapaglia@sfu.ca
% Last Update: December 2017

[m,n] = size(A);
B = zeros(m,n);
switch opt 
    case 'col'  % columns
        norms = zeros(n,1);
        for i = 1:n
            norms(i) = norm(A(:,i));
            if norms(i) ~= 0
                B(:,i) = A(:,i) / norms(i);
            else
                B(:,i) = A(:,i);
            end
        end
    case 'row'  % rows
        norms = zeros(m,1);
        for i = 1:m
            norms(i) = norms(A(i,:));
            if norms(i) ~= 0
                B(i,:) = A(i,:) / norms(i);
            else
                B(i,:) = A(i,:);
            end
        end
    otherwise
        error('parameter opt is not valid')
end
norms = norms(:);


