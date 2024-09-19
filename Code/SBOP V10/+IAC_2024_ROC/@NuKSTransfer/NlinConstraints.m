%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %%
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Radius
    r = dot(s(1:4,:), s(1:4,:), 1);
    
    % Inequality constraints
    c = [
            +dot(u(1:3,:), u(1:3,:), 1) - (params(2).^2 .* r.^8) ...     % Thrust modulation
%             -r ...
        ];
    
    % Equality constraints
    R = LegoKS.GroupAction( beta(end) );
    R = blkdiag(R,R);

    ceq = [
            reshape(params(3:10) - R * s(1:8,end), 1, []) ...
%             obj.bilinear_function(s(1:4,:), s(5:8,:)) ...
%             u(4,:) ./ r.^2 ...
          ];
end