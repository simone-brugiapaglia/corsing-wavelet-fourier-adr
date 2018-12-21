function f = lowrankfun2D(f1,f2,x,y)
% F = LOWRANKFUN2D(F1,F2,X,Y) evaluates the two-dimensional low-rank 
% function F at points X and Y. F is defined as follows:
%
%           K
%          ---
% F(x,y) = >   F1{i}(x) x F2{i}(y).
%          ---
%          i=1
%

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)
k1 = length(f1);
k2 = length(f2);
if k1 ~= k2
    error('Size of cells must agree')
end
f = 0 * f1{1}(x).* f2{1}(y);
for i = 1:k1
    f = f + f1{i}(x) .* f2{i}(y);
end
