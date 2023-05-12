%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    c = [-u(1,:).*u(2,:) s(1,:)+u(2,:) abs(u(1,:))-2 abs(u(2,:))-2];

    % Equality constraints
    ceq = [];
end