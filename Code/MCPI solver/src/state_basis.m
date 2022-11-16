%% Project: Shape-based optimization for low-thrust transfers %%
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
                B{i} = bernstein_basis(n(i),tau);
           end
        case 'Orthogonal Bernstein'
            for i = 1:length(n)
                B{i} = OB_basis(n(i),tau);
            end

        case 'Chebyshev'
            for i = 1:length(n)
                B{i} = CH_basis('first', n(i), tau);
            end

        case 'Legendre'
             for i = 1:length(n)
                B{i} = LG_basis(n(i), tau); 
             end

        case 'Hermite'
             for i = 1:length(n)
                B{i} = HT_basis(n(i), tau);
             end
             
        case 'Laguerre'
             for i = 1:length(n)
                B{i} = LR_basis(n(i), tau); 
             end

        otherwise
            error('No valid functional polynomial basis has been selected');
    end
end