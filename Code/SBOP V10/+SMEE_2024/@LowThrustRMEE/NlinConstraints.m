%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)
    % Parameters
    mu = params(1);     % Gravitational parameter

    % Timing equation 
    w = 1 + s(2,:) .* cos(t(1,:)) + s(3,:) .* sin(t(1,:));
    k = s(4,:) .* sin(t(1,:)) - s(5,:) .* cos(t(1,:));

    dtheta = sqrt(mu * s(1,:)) .* (w ./ s(1,:)).^2 + sqrt(s(1,:) / mu) .* k ./ w .* u(3,:);

    % Dimensions 
    u = u ./ dtheta;

    % Inequality constraints
    c = [
            -dtheta ...                    % The time law is monotonically increasing
            dot(u,u,1)-params(2).^2 ...    % Thrust modulation
        ];                 

    % Equalities
    ceq = [cos(tf)-cos(params(4)) sin(tf)-sin(params(4))];
end