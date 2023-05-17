%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    c = [-params(2)-u u-params(3) -s(2,:)];

    % Equality constraints
    ceq = [u-s(4,:)];
end