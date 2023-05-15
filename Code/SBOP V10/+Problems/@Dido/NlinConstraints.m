%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    c = [s([1 3],end)-params(1)];

    % Equality constraints
    ceq = [u(3,:).^2-u(2,:).^2-u(1,:).^2];
end