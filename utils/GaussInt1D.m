function I = GaussInt1D(f,h,NGpoints)
% function I = GaussInt1D(f,h)
%
% Computes the integral of f(x) over the domain [0,1] using a composite 
% Gaussian-Legendre integration formula of step h.
%
% The domain [0,1] is subdivided in a regular grid of step h. Then, on
% each subinterval of the grid, the two-point Gaussian formula is applied.
%
% The 1D formula on [-1,1] is given by
%
%   / 1         
%   |   f(x) dx ~ f(-1/sqrt(3)) + f(1/sqrt(3))  
%   /-1
%
% And the related composite 1D forumla in [0,1] is the follwing 
%
%                    N-1
%   / 1           h  ---              h                      h
%   |   f(x) dx ~ -  >   f( p_i + --------- ) + f( p_i - --------- )
%   / 0           2  ---          2 sqrt(3)              2 sqrt(3)
%                    i=0
%
% where p_i = h/2 + i*h for and N = 1/h. 
%
%
%
% Input:
%
%  f          function handle of the form @(x,y). If x,y are two vectors, 
%             f(x,y) should be a vector of the same size, containing all 
%             the evaluations 
%
%  h          step for the integration grid
%  
%  NGpoints   Number of Gaussian points to use
%
% Output:
%
%  I   Value of the integral 
%
% Simone Brugiapaglia, 2018 (simone.brugiapaglia@sfu.ca)

% choose the correct points and weights on the reference interval [-1,1]
switch NGpoints 
    case 2
        % 2 points formula
        GptsREF = [-sqrt(1/3) sqrt(1/3)];
        WeightsREF = [1 1];
    case 3
        % 3 points formula
        GptsREF = [-sqrt(3/5) 0 sqrt(3/5)];
        WeightsREF = [5/9  8/9  5/9];
    case 4
        % 4 points formula
        x1 = sqrt((3-2*sqrt(6/5))/7);
        x2 = sqrt((3+2*sqrt(6/5))/7);
        w1 = (18+sqrt(30))/36;
        w2 = (18-sqrt(30))/36;
        
        GptsREF = [-x2, -x1, x1, x2];
        WeightsREF = [w2, w1, w1, w2];
    case 5
        x1 = 0;
        x2 = 1/3*sqrt(5-2*sqrt(10/7));
        x3 = 1/3*sqrt(5+2*sqrt(10/7));
        w1 = 128/225;
        w2 = (322+13*sqrt(70))/900;
        w3 = (322-13*sqrt(70))/900;
        
        GptsREF = [-x3, -x2, x1, x2, x3];
        WeightsREF = [w3, w2, w1, w2, w3]; 
    otherwise
        error('Number of Gaussian points is not valid')
end


N = round(1/h); % number of subintervals

% create array of midpoints
midpts = h/2 : h : 1-h/2;

% create repeated array for composite integration
midptsMAT  = ones(NGpoints,1) * midpts;
midptsREP  = midptsMAT(:)';
GptsREP    = h/2 * repmat(GptsREF,1,N);

Weights1D = repmat(WeightsREF,1,N);
Gpts1D = midptsREP + GptsREP;

%plot(GaussPoints(1,:),GaussPoints(2,:),'*',[0,1,1,0],[0,0,1,1],'-')

I =  h / 2 * sum(f(Gpts1D) .* Weights1D);
end

