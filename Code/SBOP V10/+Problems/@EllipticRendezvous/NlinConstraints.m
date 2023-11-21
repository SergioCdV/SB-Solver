%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % TH transformation matri
    omega = params(4)^2 / params(6)^3;                                     % True anomaly angular velocity
    k = 1 + params(5) * cos(tau(1,:));                                     % Transformation parameter

    % Inequality constraints
    R1 = +dot(s(1:3,end), s(1:3,end))-(k(end) * params(10))^2; 
    L =  -dot(s(1:3,end), s(1:3,end))+(k(end) * params(11))^2;

    c = [
          params(9) * sqrt(dot(u(1:3,:),u(1:3,:),1))-params(3) ...         % Constraint on the force magnitude (second order cone)
          R1 ...                                                           % Keep-out sphere constraint
          L ...                                                            % Graspling reach
        ];                                                               

    % Equality constraints
    ceq = [];
end