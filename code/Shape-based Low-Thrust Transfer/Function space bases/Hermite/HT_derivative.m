%% Project: 
% Date: 30/04/22

%% Legendre derivative %%
% This function allows to compute all Hermite polynomials derivatives of order n,
% evaluated at the argument u. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated
%         - scalar degree, the degree of the derivative to be computed

% Outpus: - vector Pn, containing the evaluated Hermite polynomials
%           derivatives

function [B] = HT_derivative(order, u, degree)
    % Switch the derivative order
    switch (degree)
        case 1
            B = dhermite(order, u);
        case 2
            B = ddhermite(order, u);
        otherwise
            error('A higher-order Hermite polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Hermite tangent space
function [dPn] = dhermite(order, u)
    % Preallocation of the polynomials and its derivatives
    Pn = HT_basis(order,u);
    dPn = zeros(order+1,1); 

    % Bonnet's formula 
    for i = 2:order+1
        dPn(i) = 2*i*Pn(i-1); 
    end
end

% Second order basis of the Hermite tangent space
function [ddPn] = ddhermite(order, u)
    % Preallocation of the polynomials and its derivatives
    dPn = HT_derivative(order,u,1); 

    % Initialization of the polynomials 
    ddPn(1) = 0; 
    ddPn(2) = 0; 

    % Bonnet's formula 
    for i = 2:order
        n = i-1;
        dPn(i+1) = 2*u*ddPn(i)+4*dPn(i)-2*n*ddPn(i-1); 
    end
end