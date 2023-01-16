%% Project: 
% Date: 09/01/23

%% Legendre differentiation matrix %%
% This function allows to compute the Legendre polynomials differentiation matrix at LGL nodes. 

% Inputs: - scalar m, the number of nodes
%         - scalar u, the normalized domain value on which the polynomials
%           are to be evaluated

% Outpus: - vector D, the differentiation matrix

function [D] = LG_diffm(m, tau)
    % Compute the Legendre polynomials 
    L = LG_basis(m, tau); 

    % Preallocation 
    D = zeros(m+1); 

    % Compute the differentiation matrix
    for i = 0:m 
        for j = 0:m
            if (i ~= j)
                D(i+1,j+1) = (L(end,i+1)/L(end,j+1))/(tau(i+1)-tau(j+1));
            elseif (i == j && j == m)
                D(i+1,j+1) = (m*(m+1)/4);
            elseif (i == j && j == 0)
                D(i+1,j+1) = -(m*(m+1)/4);
            else 
                D(i+1,j+1) = 0; 
            end
        end
    end

    D = -D; 
end