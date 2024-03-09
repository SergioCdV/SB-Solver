%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    c = dot(u(1:3,:),u(1:3,:),1)-params(end)^2;         % Constraint on the force magnitude (second order cone)

    % Equality constraints
    ceq = [];
end