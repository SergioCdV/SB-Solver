%% Project: 
% Date: 01/11/2023

%% Legendre-Bernstein transformation matrix %% 
% Function to transform from the shifted Legendre basis to the Bernstein one

% Inputs:  - scalar n, the degree of the shifted Legendre series

% Outputs: - array M, the transformation matrix

function [M] = LB_tmatrix(obj, n)
    % Compute the change-of-basis matrix 
    M = zeros(n+1,n+1);

    for J = 1:n+1
        for K = 1:n+1
            total = 0;
            j = J-1;
            k = K-1;
            for i = max(0,j+k-n):min(j,k)
                total = total + (-1)^(k+i) * ( factorial(k) / (factorial(i) * factorial(k-i)) )^2 * factorial(n-k) / (factorial(n-k -(j-i)) * factorial(j-i));
            end
            M(J,K) = total / ( factorial(n) / (factorial(j) * factorial(n-j)) );
        end
    end

    M = M.';
end