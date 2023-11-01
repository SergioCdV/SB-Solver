%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %%
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Radius
    r = dot(s(1:4,:), s(1:4,:), 1);
    
    % Inequality constraints
    c = [
            +dot(u,u,1)-(params(2) .* r.^2).^2 ...     % Thrust modulation
        ];
    
    % Equality constraints
    R = [cos(beta(end)) 0 0 -sin(beta(end)); ...
        0 cos(beta(end)) sin(beta(end)) 0; ...
        0 -sin(beta(end)) cos(beta(end)) 0; ...
        sin(beta(end)) 0 0 cos(beta(end))];
    R = blkdiag(R,R);

    ceq = [
            reshape(params(3:10) - R * s(1:8,end), 1, []) ...
            u(4,:) ...
            obj.bilinear_function(s(1:4,:), s(6:9,:)) ...
          ];
end