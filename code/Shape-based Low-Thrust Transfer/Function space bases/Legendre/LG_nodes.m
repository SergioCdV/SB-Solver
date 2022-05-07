%% Project:  
% Sergio Cuevas del Valle
% Date: 20/01/20
% File: LG_nodes
% Issue: 0 
% Validated: 

%% Legendre-Gauss nodes %%
% This scripts provides the function to compute the Legendre-Gauss nodes
% for a given domain interval and polynomial degree 

% Inputs: - scalar N, the degree of the Legendre polynomial of interest
%         - scalar a, the initial domain boundary 
%         - scalar b, the final domain boundary 

% Output: - vector x, the Legendre nodes of interest

function [x] = LG_nodes(a, b, N)
    % Polynomial orders
    N = N-1;
    N(2) = N(1)+1; 
    N(3) = N(1)+2;
    
    % Initial guess using Gauss-Radau nodes
    x = linspace(-1, 1, N(2)).'; 
    y = cos((2*(0:N(1))'+1)*pi/(2*N(1)+2))+(0.27/N(2))*sin(pi*x*N(1)/N(3));
    
    % Legendre-Gauss Vandermonde Matrix
    L = zeros(N(2),N(3));
    
    % Compute the zeros of the N+1 Legendre Polynomial using the recursion relation and the Newton-Raphson method
    y0 = 2;         % Convergence value

    % Newton-Rhapson method
    while (max(abs(y-y0)) > eps)   
        L(:,1) = 1;
        L(:,2) = y;
       
        for k = 2:N(2)
            L(:,k+1) = ( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
        end
        dL = N(3)*( L(:,N(2))-y.*L(:,N(3)) )./(1-y.^2);
     
        dy = -L(:,N(3))./dL;
        y0 = y; 
        y = y0+dy;
    end
    
    % Linear map from [-1,1] to [a,b]
    y = flip(y);
    x = (a*(1-y)+b*(1+y))/2;
end