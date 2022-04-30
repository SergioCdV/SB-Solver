%% Project: 
% Date: 30/04/22

%% Chebysev basis %%
% This function allows to compute all Chebyshev polynomials derivatives of order n of both kinds,
% evaluated at the argument u. 

% Inputs: - string kind, specifying the kind of Chebyshev polynomials to
%           compute
%         - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated
%         - scalar degree, the degree of the derivative to be computed

% Outpus: - vector Pn, containing the evaluated Chebyshev polynomials 

function [B] = CH_derivative(kind, order, u, degree)
    % Switch the derivative order
    switch (degree)
        case 1
            B = dchebyshev(kind, order, u);
        case 2
            B = ddchebyshev(kind, order, u);
        otherwise
            error('A higher-order Bernstein polynomial derivative is required, but has not been implemented')
    end
end

%% Auxiliary functions 
% First order basis of the Chebyshev tangent space
function [Pn] = dchebyshev(kind, order, u)
    % Preallocation of the polynomials 
    Pn = zeros(order+1,1); 

    % Main computation 
    switch (kind)
        case 'first'
            Pn(1) = 0;                    % Initialization of the Chebyshev polynomials of the first kind
            Pn(2) = 1;                    % Initialization of the Chebyshev polynomials of the first kind

        case 'second '
            Pn(1) = 0;                    % Initialization of the Chebyshev polynomials of the second kind
            Pn(2) = 2;                    % Initialization of the Chebyshev polynomials of the second kind
        otherwise
            error('No valid kind of polynomials was selected'); 
    end
  
    for i = 2:order
        Pn(i+1) = 2*Pn(i)+2*u*Pn(i)-Pn(i-1);  % Chebyshev polynomials derivatives 
    end
end

% Second order basis of the Chebyshev tangent space
function [Pn] = ddchebyshev(kind, order, u)
    % Preallocation of the polynomials 
    Pn = zeros(order+1,1); 
    dPn = dchebyshev(kind, order, u); 

    % Main computation 
    switch (kind)
        case 'first'
            Pn(1) = 0;                    % Initialization of the Chebyshev polynomials of the first kind
            Pn(2) = 0;                    % Initialization of the Chebyshev polynomials of the first kind

        case 'second '
            Pn(1) = 0;                    % Initialization of the Chebyshev polynomials of the second kind
            Pn(2) = 0;                    % Initialization of the Chebyshev polynomials of the second kind

        otherwise
            error('No valid kind of polynomials was selected'); 
    end
  
    for i = 2:order
        Pn(i+1) = 2*Pn(i)+2*u*Pn(i)+2*dPn(i)-Pn(i-1);  % Chebyshev polynomials derivatives 
    end
end