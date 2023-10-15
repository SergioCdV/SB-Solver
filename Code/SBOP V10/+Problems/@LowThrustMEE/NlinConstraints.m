%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)
    % Constants 
    mu = params(1);                                 % Gravitational parameter of the problem
    
    % Sundman transformation
    w = 1+s(2,:).*cos(t)+s(3,:).*sin(t);
    gamma = sqrt(mu*s(1,:)).*(w./s(1,:)).^2+sqrt(s(1,:)/mu)./w.*(s(4,:).*sin(t) + s(5,:).*cos(t)) .* u(3,:);

    % Inequality constraints
    ct = dot(u,u,1)-(params(2)).^2;                 % Thrust modulation

    % Equalities (Sundman transformation)
    ceq = [cos(t(end))-cos(params(4)) sin(t(end))-sin(params(4))].';
    ceq = [];
    c = [-1./gamma.'; ct.'];
end