%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraintst(obj, params, beta, t0, tf, tau, s, u)
    % Anomaly space 
    t = (tf-t0) * 0.5 * (1+tau) + params(3);
%     t = (tf-t0) * (tau) + params(3);

    % Cross-product law 
    B = params(4) * [sin( params(2) ) * cos(t); repmat( -cos( params(2) ), 1, length(tau) );  2 * sin( params(2) ) * cos(t)];
    ub = u; 
    for i = 1:length(tau)
        ub(:,i) = cross(B(:,i), u(:,i)) / norm(B(:,i))^2;
    end

    % Inequality constraints
    sat = abs(ub) - repmat(params(14:16), 1, size(ub,2));

    c = [reshape(sat, 1, []) ...                                          % Torque saturation
         dot(s(1:3,:), s(1:3,:), 1)-params(end)^2 ...                     % Angular velocity tolerance
         ]; 

    for i = 1:length(tau)
        ub(:,i) = cross(ub(:,i), B(:,i));
    end

    c = [c dot(ub-u, ub-u, 1) - 1e-8];

    % Equality constraints
    ceq = [];
end