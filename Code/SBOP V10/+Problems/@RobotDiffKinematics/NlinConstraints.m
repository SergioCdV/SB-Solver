%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    c = [-reshape(u(1:3,:), 1, [])-params(3) ... % Constraint on the torque magnitude (infinity norm in epigraph form)
          reshape(u(1:3,:), 1, [])-params(3) ... % Constraint on the torque magnitude (infinity norm in epigraph form)
         -reshape(u(4:6,:), 1, [])-params(4) ... % Constraint on the torque magnitude (infinity norm in epigraph form)
          reshape(u(4:6,:), 1, [])-params(4) ... % Constraint on the torque magnitude (infinity norm in epigraph form)
        ];                                                                                      

    % Equality constraints
    ceq = [];    
end