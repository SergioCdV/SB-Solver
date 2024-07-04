%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)
    % Parameters
    mu = params(1);     % Gravitational parameter

    % Timing equation 
    l = t(1,:);
    w = 1 + s(2,:) .* cos(l) + s(3,:) .* sin(l);
    k = s(4,:) .* sin(l) - s(5,:) .* cos(l);

    dtheta = sqrt(mu * s(1,:)) .* (w ./ s(1,:)).^2 + sqrt(s(1,:) / mu) .* k ./ w .* u(3,:);

    % Inequality constraints
    c = [
            -dtheta ...                     % The time law is monotonically increasing
            dot(u,u,1) - params(2)^2 ...    % Thrust modulation
        ];                 

    % Equalities
%     ceq = [cos(t(1,end))-cos(params(4)) sin(t(1,end))-sin(params(4))];
%     ceq = [dot(t(2,:), 1./dtheta) - params(6)];
    ceq = [];
end