%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %%
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Radius
%     mu = params(1); 
%     [~, alpha] = LegoKS.OscEnergy(mu, s, "Ecc");
    r = dot(s(1:4,:), s(1:4,:), 1);
    
    % Inequality constraints
    c = [
%             +dot(u(1:3,:), u(1:3,:), 1) - (params(2).^2 .* r.^4) ...     % Thrust modulation
            r(1) - 1.001 * r ...                                         % Bound on the radial distance
            r - 1.001 * r(end) ...                                       % Bound on the radial distance
        ];
    
    % Equality constraints
    ceq = [
            cos(tau(1,end)) - cos(params(end-1)) ...
            sin(tau(1,end)) - sin(params(end-1)) ...
          ];
end