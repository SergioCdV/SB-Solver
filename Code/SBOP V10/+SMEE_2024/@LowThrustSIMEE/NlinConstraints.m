%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)
    % Parameters
%     mu = params(1);     % Gravitational parameter

    % Compute the evolution of the equinoctial latitude
%     l = t(1,:) + s(6,:);
%     w = 1 + s(2,:) .* cos(l) + s(3,:) .* sin(l);
%     dtheta = sqrt(mu * s(1,:)) .* (w ./ s(1,:)).^2; 

    % Inequality constraints
    u_square = dot( u(1:3,:), u(1:3,:), 1 );
    c = [
             u_square - params(2)^2
%             +u_square - 1.01 ...               % Thrust modulation
%             -u_square + 0.99 ...               % Thrust modulation
%             +dot(t(2,:), 1 ./ dtheta) - 1.01 * params(6) ...
%             -dot(t(2,:), 1 ./ dtheta) + 0.99 * params(6) ...
        ];                 

    % Equalities
    ceq = [cos(t(1,end) + beta)-cos(params(4)) sin(t(1,end) + beta)-sin(params(4))];
end