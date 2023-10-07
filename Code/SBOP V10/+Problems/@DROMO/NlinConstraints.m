%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)    
    % Sundman transformation
    S = 1+s(2,:).*cos(t)+s(3,:).*sin(t);
    gamma = s(3,:).^3.*S.^2;

    % Inequality constraints
    ct = dot(u,u,1)-(params(2) * gamma.*s(3,:).*S).^2;   % Thrust modulation

    % Equalities (Sundman transformation)
    % ceq = [cos(t(end))-cos(params(3)) sin(t(end))-sin(params(3))].';
    % ceq = [];
    c = [-1./gamma.'; ct.'];
end