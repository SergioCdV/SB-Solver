%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    c = [s(1,:)-params(2) s(2,:)-params(3) s(4,:)-params(4) dot(u,u,1)-params(5)^2];
   
    % Equality constraints
    ceq = [];
end