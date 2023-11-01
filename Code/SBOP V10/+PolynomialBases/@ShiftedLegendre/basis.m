%% Project: 
% Date: 01/11/23

%% Shifted Legendre basis %%
% This function allows to compute all shifted Legendre polynomials of order n,
% evaluated at the argument u. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector Pn, containing the evaluated Legendre polynomials 

function [Pn] = basis(obj, order, u)
    % Preallocation of the polynomials 
    Pn = zeros(order+1,length(u)); 
    u = 2 * u - 1;

    % Initialization of the polynomials 
    Pn(1,:) = ones(1,length(u)); 
    Pn(2,:) = u; 

    % Bonnet's formula 
    for i = 3:order+1
        n = i-1;
        Pn(i,:) = ((2*n-1)*u.*Pn(i-1,:)-(n-1)*Pn(i-2,:))/n; 
    end
end
