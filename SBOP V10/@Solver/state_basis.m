%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/04/22

%% State basis %%
% Function to preallocate the state vector polynomial use of interest

% Inputs: - scalar L, the degree of the greatest derivative in the dynamics
%         - scalar n, the maximum polynomial order of the approximation
%         - string basis, the polynomial basis to be used
%         - vector tau, the control parameter vector 

% Outputs: - cell array B, the state vector shape base basis

function [B, tau] = state_basis(obj, L, n, basis, tau)
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
        B{i} = P.basis(n(i),tau);
        for j = 1:L
            B{i} = [B{i}; P.derivative(n(i), tau, j)];
        end
    end
end