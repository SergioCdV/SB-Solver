%% Project: 
% Date: 30/01/2022

%% Bernstein derivative %% 
% Function to compute the basis of the tangent space of the Bernstein polynomials 

% Inputs:  - scalar n, the cardinal of the polynomial basis
%          - vector tau, the parametrization of the basis
%          - scalar order, the order of the derivative required

% Outputs: - array B, the basis of the required tangent space, of size n+1 x length(tau)

function [B] = bernstein_derivative(n, tau, order)
    % Switch the derivative order
    switch (order)
        case 1
            B = dbernstein(n, tau);
        case 2
            B = ddbernstein(n, tau);
        otherwise
            error('A higher-order Bernstein polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Bernstein tangent space
function [B] = dbernstein(n, tau)
    % Preallocation for speed 
    B = zeros(n+1,length(tau));

    % Computation of the basis 
    for i = 0:n
        % Indexing parameter 
        j = i+1; 

        % Definition of the basis
        if (i == 0)
            B(j,:) = -n.*(1-tau).^(n-1);
            
        elseif (i < n)
            K(1) = factorial(n) / (factorial(i-1)*factorial(n-i));
            K(2) = factorial(n) / (factorial(i)*factorial(n-i-1));

            B(j,:) = K(1)*tau.^(i-1).*(1-tau).^(n-i) - K(2)*tau.^i.*(1-tau).^(n-i-1);

        else
            B(j,:) = n.*tau.^(n-1);
        end
    end
end

% Second order basis of the Bernstein tangent space 
function [B] = ddbernstein(n, tau)
    % Preallocation for speed 
    B = zeros(n+1,length(tau));

    % Computation of the basis 
    for i = 0:n
        % Indexing parameter 
        j = i+1; 

        % Definition of the basis
        if (i == 0)
            B(j,:) = n.*(n-1).*(1-tau).^(n-2);

        elseif (i == 1)
            B(j,:) = n.*(n-1).*(n-2).*tau.*(1-tau).^(n-3)-2.*n.*(n-1).*(1-tau).^(n-2);
            
        elseif (i == (n-1))
            B(j,:) = n.*(n-1).*(n-2).*tau.^(n-3).*(1-tau)-2.*n.*(n-1).*tau.^(n-2);

        elseif (i == n)
            B(j,:) = n.*(n-1).*tau.^(n-2);
           
        else
            K(1) = factorial(n) / (factorial(i-2)*factorial(n-i));
            K(2) = 2 * factorial(n) / (factorial(i-1)*factorial(n-i-1));
            K(3) = factorial(n) / (factorial(i)*factorial(n-i-2));

            B(j,:) = K(1)*tau.^(i-2).*(1-tau).^(n-i) - K(2)*tau.^(i-1).*(1-tau).^(n-i-1) + K(3)*tau.^i.*(1-tau).^(n-i-2);
        end
    end
end