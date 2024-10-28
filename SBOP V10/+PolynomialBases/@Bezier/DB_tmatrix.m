%% Project: 
% Date: 01/11/2023

%% Derivative Bernstein transformation matrix %% 
% Function to transform from the derivative of the Berstein polynomial coefficients

% Inputs:  - scalar n, the degree of the shifted Legendre series

% Outputs: - array M, the transformation matrix

function [M] = DB_tmatrix(obj, n, order)
    % Compute the change-of-basis matrix 
    M = zeros(n+1,n+1);

%     n = 3;
%     q = 1;
% 
%     [a_0 a_1 a_2 a_3] = [a_0 a_1 a_2 a_3] * [C-1(0,n,q) C0(1,n,q) 0 * C1(2,n,q) 0 * C(3,n,q); ...
%                                              C-0(0,n,q) C1(1,n,q) C2(2,n,q) C2(3,n,q); ...
%                                              0 * C+1(0,n,q) C2(1,n,q) C3(2,n,q) C3(3,n,q) ...
%                                              0 * (0,n,q) 0 * C(1,n,q) C(2,n,q) C4(3,n,q)];

%     n = 2;
%     q = 1;
% 
%     [a_0 a_1 a_2] = [a_0 a_1 a_2] * [         0 C0(1,n,q) C1(2,n,q)
%                                      C-0(0,n,q) C1(1,n,q) C2(2,n,q)
%                                      C+1(0,n,q) C2(1,n,q) 0];

    for J = 1:n+1
        j = J-1;
        K = max(0,j-order):min(n,j+order);
        total = zeros(length(K),1);

        for k = 1:length(K)
            q = K(k);
            for i = 0:order
                if (j - i - q >= 0 && n - j - order + i + q >= 0 && order - i - q >= 0)
                    C = (-1)^(order+i) * ( factorial(order) / (factorial(i) * factorial(order - i)) );
                    C = C * ( factorial(j) / (factorial(q + i) * factorial(j - i - q)) );
                    total(k) = total(k) + C * ( factorial(n - j) / (factorial(order - i - q) * factorial(n - j - order + i + q)) );
                end
            end
        end

        M(K+1,J) = total;
       
    end

    M = factorial(order) * M;
end