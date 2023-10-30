%% Project: 
% Date: 30/04/22

%% Legendre modal projection %%
% This function allows to project the trajectory onto the basis of Legendre polynomials. 

% Inputs: - scalar order, determining the order of the approximation 
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated
%         - scalar degree, the degree of the derivative to be computed

% Outpus: - vector Pn, containing the evaluated Legendre polynomials
%           derivatives

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

