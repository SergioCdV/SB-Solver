%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Mayer and Lagrange terms
    if (length(params) > 58)
        M = dot( s(1:obj.StateDim,end) - params(59:64).', s(1:obj.StateDim,end) - params(59:64).' ); 
    else
        M = 0;
    end

    L = dot(s(7:12,:), s(7:12,:), 1);
end