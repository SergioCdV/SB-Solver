%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)
    % Constants 
    mu = params(1);                                 % Gravitational parameter of the problem
    
    % Sundman transformation
    l = s(6,:);
    w = 1 + s(2,:) .* cos(l) + s(3,:) .* sin(l);
    k = s(4,:) .* sin(l) - s(5,:) .* cos(l);

    dtheta = sqrt(mu * s(1,:)) .* (w ./ s(1,:)).^2 + sqrt(s(1,:) / mu) .* k ./ w .* u(3,:);

    % Inequality constraints
    u_square = dot( u(1:2,:), u(1:2,:), 1) / params(2)^2;
    c = [
            -dtheta ...
            +u_square - 1.02 ...               % Thrust modulation
            -u_square + 0.98 ...               % Thrust modulation
%             +dot(t(2,:), 1 ./ dtheta) - 1.01 * params(6) ...
%             -dot(t(2,:), 1 ./ dtheta) + 0.99 * params(6) ...
        ];


    % Equalities (Sundman transformation)
    ceq = [cos(t(end))-cos(params(4)) sin(t(end))-sin(params(4))].';
    ceq = [s(12,:) - dtheta];
end