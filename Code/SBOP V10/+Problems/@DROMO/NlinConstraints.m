%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)    
    % Sundman transformation
    S = 1+s(1,:).*cos(t)+s(2,:).*sin(t);
    gamma = s(3,:).^3.*S.^2;

    % Inequality constraints
    ct = dot(u,u,1)-params(2).^2;   % Thrust modulation

    % Equalities (Sundman transformation)
    ceq = [cos(t(end))-cos(params(4)) sin(t(end))-sin(params(4)) dot(s(4:7,:),s(4:7,:),1)-1].';
    c = [-gamma ct -s(3,:)].';
end