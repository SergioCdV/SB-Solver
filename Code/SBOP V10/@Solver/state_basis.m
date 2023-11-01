%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/04/22

%% State basis %%
% Function to preallocate the state vector polynomial use of interest

% Inputs: - scalar L, the degree of the greatest derivative in the dynamics
%         - scalar n, the maximum polynomial order of the approximation
%         - string basis, the polynomial basis to be used
%         - vector tau, the control parameter vector 

% Outputs: - cell array B, the state vector shape base basis
%          - cell array C, the state vector shape base basis to be used to
%            evaluate constraints

function [B, C, tau] = state_basis(obj, L, n, basis, tau)
    % Extract the polynomial order 
    if (~exist('n', 'var'))
        n = obj.PolOrder;
    end
    
    if (~exist('n', 'var'))
        basis = obj.Basis;
    end
    
    if (~exist('tau', 'var'))
        tau = obj.gridding().tau;
    end

    % Preallocation of the polynomial state vector support
    B = cell(length(n),1);   
    C = cell(length(n),1); 

    switch (basis)
        case 'Legendre'
            ctau = 0.5 * (tau + 1);
        otherwise
            ctau = tau;
    end

    switch (basis)
        % Finite-horizon polynomial basis
        case 'Bernstein'
            P = PolynomialBases.Bezier();

        case 'Orthogonal Bernstein'
            P = PolynomialBases.OBezier();

        case 'Chebyshev'
            P = PolynomialBases.Chebyshev();

        case 'Legendre'
            P = PolynomialBases.Legendre();

        % Infinite-horizon polynomial basis (not currently supported)
        case 'Hermite'
            P = PolynomialBases.Hermite();
             
        case 'Laguerre'
            P = PolynomialBases.Laguerre();

        otherwise
            error('No valid functional polynomial basis has been selected');
    end

    for i = 1:length(n)
        base = P.basis(n(i),tau);

        derivative = zeros(L * (n(i)+1),length(tau));
        cderivative = zeros(L * (n(i)+1),length(ctau));
        
        switch (basis)
            case 'Legendre'
                M = PolynomialBases.Bezier().LB_tmatrix(n(i));
                cbase = M.' * PolynomialBases.Bezier().basis(n(i), ctau);
            otherwise
                cbase = base;
        end

        for j = 1:L
            derivative(1+(n(i)+1)*(j-1):(n(i)+1)*j,:) = P.derivative(n(i), tau, j);

            switch (basis)
                case 'Legendre'
                    cderivative(1+(n(i)+1)*(j-1):(n(i)+1)*j,:) = M.' * PolynomialBases.Bezier().derivative(n(i), ctau, j);
                otherwise
                    cderivative(1+(n(i)+1)*(j-1):(n(i)+1)*j,:) = derivative(1+(n(i)+1)*(j-1):(n(i)+1)*j,:);
            end
        end

        B{i} = [base; derivative];
%         C{i} = [cbase; cderivative];
        C{i} = B{i};
 
    end
end