%% Project: 
% Date: 30/10/23

%% Chebyshev modal projection %%
% This function allows to project the trajectory onto the basis of Chebyshev polynomials. 

% Inputs: - array S, the trajectory to be projected

% Outpus: - array C, the array of modal coefficients

function [C] = modal_projection(obj, S)
    % Constants 
    n = size(S,2);      % Dimension of the projection

    % Compute the nodes and the polynomials
    tau = cos(pi * ((0:n-1)+0.5) / n);
    P = obj.basis(n-1, tau);
    
    % Compute the coefficients 
    C = zeros(size(S));

    for i = 1:n
        C(:,i) = (2 - (i == n)) / n .* sum(P(i,:) .* S, 2); 
    end
end

