function answer = isconst(f,dom)
%ANSWER = ISCONST(F,DOM) checks whether a function handle F defined over a 
% domain DOM to C is constant or not.
%
% INPUT:
%   F      function handle F : [a,b] --> C
%   DOM    interval [a,b] (default = [0,1])
%
% OUTPUT:
%   ANSWER  1 if f is constant, 0 otherwise 
%

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

if nargin == 1
    dom = [0,1];
end

% domain
a = dom(1); 
b = dom(2);

vals   = a + (b-a) * rand(1000,1); % generate random values
g      = @(x) x*0 + f(vals(1));    % constant function f(val(1))
answer = prod(f(vals) == g(vals)); % f == f(val(1)) ?

end

