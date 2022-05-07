%% Project:  
% Sergio Cuevas del Valle
% Date: 20/01/20
% File: LG_nodes
% Issue: 0 
% Validated: 

%% Laguerre nodes %%
% This scripts provides the function to compute the Laguerre nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Legendre polynomial of interest
%         - scalar a, the initial domain boundary 
%         - scalar b, the final domain boundary 

% Output: - vector x, the Laguerre nodes of interest

function [x] = LR_nodes(a, b, N)
    % Polynomial orders
    N = N-1;
    N(2) = N(1)+1; 
    N(3) = N(2)+1;
    
    % Initial guess using Gauss-Radau nodes
    x = linspace(0, 1e3, N(2)).'; 
    y = cos((2*(0:N(1))'+1)*pi/(2*N(1)+2))+(0.27/N(2))*sin(pi*x*N(1)/N(3));
    
    % Legendre-Gauss Vandermonde Matrix
    L = zeros(N(2),N(3));
    dL = zeros(N(2),N(3));

    % Compute the zeros of the N+1 Laguerre Polynomial using the recursion relation and the Newton-Raphson method
    y0 = 2;         % Convergence value

    % Index 
    while (max(abs(y-y0)) > eps)   
        for i = 1:length(y)
            L(i,:) = LR_basis(N(2),y(i)).';
            dL(i,:) = LR_derivative(N(2),y(i),1).';
        end

        dy = -L(:,N(3))./dL(:,N(3));
        y0 = y;
        y = y0+dy;
        max(abs(y-y0))
    end
    
    % Linear map from [-1,1] to [a,b]
    x = (a*(1-y)+b*(1+y))/2;
end