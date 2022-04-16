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
    N(2) = N(1)+1; 
    N(3) = N(2)+1;
    
    % Initial guess using Gauss-Randau nodes
    y = -cos(2*pi*(0:N(1))/(2*N(1)+1))';
    
    % Legendre-Gauss Vandermonde Matrix
    L = zeros(N(2),N(3));
    
    % Compute the zeros of the N+1 Legendre Polynomial using the recursion relation and the Newton-Raphson method
    y0 = 2;         % Convergence value

    % Index 
    index = 2:N(2); 
    while (max(abs(y-y0)) > eps)   
        y0 = y;
        L(1,:)=(-1).^(0:N(2));       
        L(index,1) = 1;    
        L(index,2) = y(index);
       
        for k = 2:N(2)
            L(index,k+1) = ( (2*k-1)*y(index).*L(index,k)-(k-1)*L(index,k-1) )/k;
        end
     
        dy = ((1-y0(index))/N(2)) .* (L(index,N(2))+L(index,N(3))) ./ (L(index,N(2))-L(index,N(3)) );
        y(index) = y0(index)-dy;
    end
    
    % Linear map from [-1,1] to [a,b]
    x = (a*(1-y)+b*(1+y))/2;
end