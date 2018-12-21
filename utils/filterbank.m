function [primal,dual] = filterbank(wname) 
% [PRIMAL,DUAL] = FILTERBANK(WNAME) outputs two structs PRIMAL and DUAL 
% containing the primal and dual wavelet filters associated with the 
% wavelet basis WNAME.
%
% WNAME is a string and it can be 
%  'bior' biorthogonal B-splines of order d = 2 and dtilde = 2,
%  'hier' hiearchical B-splines of order d = 2.
%
% The structs are as follows:
%
% PRIMAL.a
% PRIMAL.b
% PRIMAL.a_start
% PRIMAL.b_start
%
% DUAL.a
% DUAL.b
% DUAL.a_start
% DUAL.b_start
%

% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)


switch wname
    case 'bior'
        primal.a = [0.5 1 0.5];
        primal.a_start = -1;
        
        primal.b = [0.25 0.5 -1.5 0.5 0.25];
        primal.b_start = -1;
        
        dual.a = [-1/4, 1/2, 3/2, 1/2, -1/4];
        dual.a_start = -2;
        
        dual.b = [1/2, -1, 1/2];
        dual.b_start = 0;
        
    case 'hier'
        primal.a = [0.5 1 0.5];
        primal.a_start = -1;
        
        primal.b = [1];
        primal.b_start = 1;
        
        dual.a = [];
        dual.a_start = -Inf; 
        
        dual.b = [];
        dual.b_start = -Inf; 
        
    otherwise
        error('Trial family is not valid')
end

end

