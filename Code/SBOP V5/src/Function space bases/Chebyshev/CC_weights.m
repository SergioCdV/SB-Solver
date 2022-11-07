%% Project:  
% Sergio Cuevas del Valle
% Date: 06/11/22
% File: CC_weights
% Issue: 0 
% Validated: 

%% Clenshaw-Curtis weights %%
% This scripts provides the function to compute the Clenshaw-Curtis weights
% for a given domain interval and polynomial degree 

% Inputs: - vector x, the Legendre nodes at which the weights shall be
%           evaluated
%         - vector dP, the n-th Legendre derivative at the Legendre nodes

% Output: - vector w, the Legendre weights of interest

function [w] = CC_weights(N)  
    % DCT-I matrix 
    D = zeros(N/2+1);                           % Preallocation 

    for i = 1:(N/2+1)
        for j = 1:(N/2+1)
            if ((j-1) == 0 || ((j-1) == N/2))
                D(i,j) = cos((j-1)*(i-1)*pi/(N/2))/N;
            else
                D(i,j) = 2*cos((j-1)*(i-1)*pi/(N/2))/N;
            end
        end
    end

    % Final computation
    d = [1; 2./(1-(2:2:N-2).'.^2); 1/(1-N^2)];      % Coefficients
    w = D.'*d;                                      % Final weights
end