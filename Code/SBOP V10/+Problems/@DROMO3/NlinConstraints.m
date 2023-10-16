%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u) 
    % Constants 
    l = 2*s(1,:)-s(1,:).^2-s(2,:).^2;

    % Inequality constraints
    ct = dot(u,u,1)-params(2).^2;   % Thrust modulation
    c = [-s(3,:).*l ct].';

    % Equalities 

    ceq = [dot(s(4:7,:),s(4:7,:),1)-1 ...                       % Quaternion norm constraint
           l(end)-(1-params(3)^2) ...                           % Final eccentricity
           s(12,:)-2*s(3,:).^3.*s(1,:).*l.*(u(1,:)+u(2,:));     % Dynamic constraint
           ].';
end