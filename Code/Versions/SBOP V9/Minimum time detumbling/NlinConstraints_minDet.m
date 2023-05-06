%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    c = [dot(u,u,1)-params(1)^2 abs(s(4,:))-params(2) abs(s(5,:))-params(2) abs(s(6,:))-params(2)]; 

    % Equality constraints
    ceq = [];
end