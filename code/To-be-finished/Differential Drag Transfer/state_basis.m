%% Project: 
% Date: 16/04/22

%% State basis %%
% Function to preallocate the state vector polynomial use of interest

% Inputs: - scalar n, the degree of the approximation 
%         - vector tau, the control parameter vector 
%         - string basis, the Bernstein polynomial basis to be used

% Outputs: - cell array B, the state vector shape base basis

function [B] = state_basis(n, tau, basis)
    % Preallocation of the Bernstein basis
    B = cell(length(n),1);              

    % Reshape the Bernstein basis
    for i = 1:length(n)
        switch (basis)
            case 'Bernstein'
                B{i} = [bernstein_basis(n(i),tau); bernstein_derivative(n(i),tau,1); bernstein_derivative(n(i),tau,2)];
            case 'Orthogonal Bernstein'
                B{i} = [OB_basis(n(i),tau); OB_derivative(n(i),tau,1); OB_derivative(n(i),tau,2)];
            otherwise
                error('No valid collocation polynomial basis has been selected');
        end
    end
end