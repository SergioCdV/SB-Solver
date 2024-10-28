%% Project: 
% Date: 30/10/23

%% Legendre modal projection %%
% This function allows to project the trajectory onto the basis of Legendre polynomials. 

% Inputs: - array S, the trajectory to be projected

% Outpus: - array C, the array of modal coefficients

function [C] = modal_projection(obj, S)
    % Constants 
    n = size(S,2);      % Dimension of the projection

    % Compute the nodes and the polynomials
    grid = CollocationMesh.LegendreGrid(n-1);
    P = obj.basis(n-1, grid.tau);
    
    % Compute the coefficients 
    C = zeros(size(S));

    for i = 1:n
        C(:,i) = ( i + 0.5) .* sum(grid.W .* P(i,:) .* S, 2); 
    end
end

