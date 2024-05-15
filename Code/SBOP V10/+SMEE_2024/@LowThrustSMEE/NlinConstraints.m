%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)
    % Parameters
    mu = params(1);     % Gravitational parameter

    % Compute the evolution of the equinoctial latitude
    l = t(1,:) + s(6,:);
    w = 1 + s(2,:) .* cos(l) + s(3,:) .* sin(l);
    dtheta = sqrt(mu * s(1,:)) .* (w ./ s(1,:)).^2; 

    % Inequality constraints
    c = [
            -dtheta; ...                                % Sundman transformation positiveness
            dot(u,u,1)-(params(2)).^2 ...               % Thrust modulation
        ];                 

    % Equalities
    ceq = [cos(l(end))-cos(params(4)) sin(l(end))-sin(params(4))];
end