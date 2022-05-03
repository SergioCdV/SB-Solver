%% Project: 
% Date: 30/04/22

%% Legendre basis %%
% This function allows to compute all Legendre polynomials of order n of both kinds,
% evaluated at the argument u. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector Pn, containing the evaluated Chebyshev polynomials 

function [Pn] = LG_basis(order, u)
    % Preallocation of the polynomials 
    Pn = zeros(order+1,1); 

    for i = 1:order+1
        n = i-1;
        for k = 0:n
            Pn(i) = Pn(i)+(1/2^n)*(factorial(n)/(factorial(k)*factorial(n-k)))^2*(u+1)^(n-k)*(u-1)^k; 
        end
    end
end
