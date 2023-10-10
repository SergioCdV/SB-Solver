%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % TH transformation matrix
    omega = params(4)^2 / params(6)^3;                                     % True anomaly angular velocity
    k = 1 + params(5) * cos(tau(1,:));                                     % Transformation parameter
    kp =  - params(5) * sin(tau(1,end));                                   % Derivative of the transformation parameter
    L = [k(end) * eye(3) zeros(3); kp * eye(3) eye(3)/(k(end) * omega)];   % TH transformation matrix
    SF = L \ s(1:6,end);                                                   % Physical boundary conditions

    % Inequality constraints
    R = dot(params(9:11).'-SF(1:3)-params(12:14).', params(9:11).'-SF(1:3)-params(12:14).', 1)-params(2)^2; 

    c = [dot(u(1:3,:),u(1:3,:),1)-params(3)^2 ...                          % Constraint on the force magnitude (second order cone)
         -k.^2 * omega ...                                                 % Monotony of the time law
         R];                                                               % Graspling constraint

    % Equality constraints
    ceq = [cos(tau(1,end))-cos(params(7)); ...                             % Multi-revolution anomaly constraint
           sin(tau(1,end))-sin(params(7)); ...                             % Multi-revolution anomaly constraint
          ];
end