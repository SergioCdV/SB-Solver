%% Project:  
% Sergio Cuevas del Valle
% Date: 07/05/20
% File: CH_nodes
% Issue: 0 
% Validated: 

%% Chebyshev nodes %%
% This scripts provides the function to compute the Chebysehv nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar a, the initial domain boundary 
%         - scalar b, the final domain boundary 
%         - scalar N, the degree of the Chebyshev polynomial of interest

% Output: - vector x, the Chebyshev nodes of interest

function [x] = CH_nodes(a, b, N)
    % Chebyshev nodes 
    i = N:-1:1;
    y = cos((2*i-1)/(2*N)*pi);
    
    % Linear map from [-1,1] to [a,b]
    x = (a*(1-y)+b*(1+y))/2;
end