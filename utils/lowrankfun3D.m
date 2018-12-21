function f = lowrankfun3D(f1,f2,f3,x,y,z)
% F = LOWRANKFUN2D(F1,F2,F3,X,Y,Z) evaluates the three-dimensional low-rank
% function F at points X, Y, and Z. F is defined as follows:
%
%             K
%            ---
% F(x,y,z) = >   F1{i}(x) x F2{i}(y) x F3{i}(z).
%            ---
%            i=1
%

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)
k1 = length(f1);
k2 = length(f2);
k3 = length(f3);
if k1 ~= k2 || k1 ~= k3
    error('Length of cells must agree')
end
f = 0 * f1{1}(x).* f2{1}(y) .* f3{1}(z);
for i = 1:k1
    f = f + f1{i}(x) .* f2{i}(y) .* f3{i}(z);
end
