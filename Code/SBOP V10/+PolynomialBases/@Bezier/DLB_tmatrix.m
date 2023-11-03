%% Project: 
% Date: 01/11/2023

%% Derivative Legendre-Bernstein transformation matrix %% 
% Function to transform from the derivative of the shifted Legendre basis to the Bernstein one

% Inputs:  - scalar n, the degree of the shifted Legendre series

% Outputs: - array M, the transformation matrix

function [M] = DLB_tmatrix(obj, n)
    % Compute the change-of-basis matrix 
    M1 = zeros(n+1,n+1);
    M2 = zeros(n+1,n+1);

    for J = 1:n+1
        for K = 1:n+1
            total = 0;
            j = J-1;
            k = K-1;
            for i = max(1,j+1+k-n):min(j+1,k)
                C = (-1)^(k+i) * k * ( factorial(k) / (factorial(i) * factorial(k-i)) );
                total = total + C * ( factorial(k-1) / (factorial(i-1) * factorial(k-i)) ) * factorial(n-k+1) / (factorial(n-k+1-(j-i+1)) * factorial(j-i+1));
            end
            M1(J,K) = total / ( factorial(n) / (factorial(j) * factorial(n-j)) );
        end
    end

    for J = 1:n+1
        for K = 1:n+1
            total = 0;
            j = J-1;
            k = K-1;
            for i = max(0,j+k-1-n):min(j,k-1)
                C = (-1)^(k+i) * k * ( factorial(k) / (factorial(i) * factorial(k-i)) );
                total = total + C * ( factorial(k-1) / (factorial(i) * factorial(k-1-i)) ) * factorial(n-k+1) / (factorial(n-k+1-(j-i)) * factorial(j-i));
            end
            M2(J,K) = total / ( factorial(n) / (factorial(j) * factorial(n-j)) );
        end
    end

    M = M1 - M2;
    M = M.';
end