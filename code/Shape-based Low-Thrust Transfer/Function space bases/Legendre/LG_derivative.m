%% Project: 
% Date: 30/04/22

%% Legendre derivative %%
% This function allows to compute all Legendre polynomials derivatives of order n of both kinds,
% evaluated at the argument u. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated
%         - scalar degree, the degree of the derivative to be computed

% Outpus: - vector Pn, containing the evaluated Chebyshev polynomials 

function [B] = LG_derivative(order, u, degree)
    % Switch the derivative order
    switch (degree)
        case 1
            B = dlegendre(order, u);
        case 2
            B = ddlegendre(order, u);
        otherwise
            error('A higher-order Legendre polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Chebyshev tangent space
function [dPn] = dlegendre(order, u)
    % Preallocation of the polynomials and its derivatives
    Pn = LG_basis(order,u);
    dPn = zeros(order+1,1); 

    % Main computation
    for i = 1:order
        dPn(i+1) = (i+1)*Pn(i)+u*dPn(i);  % Legendre polynomials derivatives 
    end
end

% Second order basis of the Chebyshev tangent space
function [ddPn] = ddlegendre(order, u)
    % Preallocation of the polynomials and its derivatives
    dPn = LG_derivative(order,u,1); 
    ddPn = zeros(order+1,1);

    % Main computation
    for i = 1:order
        ddPn(i+1) = (i+2)*dPn(i)+u*ddPn(i);  % Legendre polynomials derivatives 
    end
end