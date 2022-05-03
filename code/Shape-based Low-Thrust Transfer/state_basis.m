%% Project: 
% Date: 16/04/22

%% State basis %%
% Function to preallocate the state vector polynomial use of interest

% Inputs: - scalar n, the degree of the approximation 
%         - vector tau, the control parameter vector 
%         - string basis, the Bernstein polynomial basis to be used

% Outputs: - cell array B, the state vector shape base basis

function [B, tau] = state_basis(n, tau, basis)
    % Preallocation of the Bernstein basis
    B = cell(length(n),1);              

    % Reshape the Bernstein basis
    switch (basis)
        case 'Bernstein'
           for i = 1:length(n)
                B{i} = [bernstein_basis(n(i),tau); bernstein_derivative(n(i),tau,1); bernstein_derivative(n(i),tau,2)];
           end
        case 'Orthogonal Bernstein'
            for i = 1:length(n)
                B{i} = [OB_basis(n(i),tau); OB_derivative(n(i),tau,1); OB_derivative(n(i),tau,2)];
            end

        case 'Chebyshev'
            tau = 2*tau-1;
            for i = 1:length(n)
                for j = 1:length(tau)
                    B{i}(:,j) = [CH_basis('first', n(i), tau(j)); CH_derivative('first', n(i), tau(j), 1); CH_derivative('first', n(i), tau(j), 2)];
                end
            end

        case 'Legendre'
             tau = 2*tau-1;
             for i = 1:length(n)
                for j = 1:length(tau)
                    B{i}(:,j) = [LG_basis(n(i), tau(j)); LG_derivative(n(i), tau(j), 1); LG_derivative(n(i), tau(j), 2)];
                end
            end
        otherwise
            error('No valid collocation polynomial basis has been selected');
    end
end