%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % TH transformation matri
    omega = params(2)^2 / params(4)^3;        % True anomaly angular velocity
    k = 1 + params(3) * cos(tau(1,:));        % Transformation parameter

    % Inequality constraints
    c = dot(u(1:3,:), u(1:3,:), 1) - params(7)^2;         % Constraint on the force magnitude (second order cone)

    % Equality constraints
    ceq = [];
end