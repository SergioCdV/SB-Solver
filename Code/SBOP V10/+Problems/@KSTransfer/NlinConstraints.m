%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Radius 
    r = dot(s(1:4,:),s(1:4,:),1); 

    U = s(1,:) .* u(1,:) +s(2,:) .* u(2,:) + s(3,:) .* u(3,:) + s(4,:) .* u(4,:); 

    % Inequality constraints
    c = [dot(u,u,1)-params(2)^2 .* r ...
        -dot(s(1:4,:),s(1:4,:),1)];

    % Equality constraints
    ceq = [s(10,:)+4/params(1).*dot(s(6:9,:),r.*u,1) ... 
           obj.bilinear_function(s(1:4,:), s(6:9,:)) ...
           U];
end