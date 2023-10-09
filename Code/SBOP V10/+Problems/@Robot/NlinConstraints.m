%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % TH transformation matrix
    omega = params(23)^2 / params(25)^3;                                    % True orbit mean motion
    k = 1 + params(24) * cos(tau(1,:));                                     % Transformation parameter
    kp =  - params(24) * sin(tau(1,end));                                   % Derivative of the transformation parameter
    L = [k(end) * eye(3) zeros(3); kp * eye(3) eye(3)/(k(end) * omega)];    % TH transformation matrix
    SF = L \ s([1:3 8:10],end);                                             % Physical boundary conditions

    % Inequality constraints
    R = dot(params(4:6).'-SF(1:3), params(4:6).'-SF(1:3), 1)-params(3)^2; 

    c = [dot(u(1:3,:),u(1:3,:),1)-params(2)^2 ...                           % Constraint on the force magnitude (second order cone)
         dot(u(4:6,:),u(4:6,:),1)-params(22)^2 ...                          % Constraint on the torque magnitude (second order cone)
         dot(s(4:7,:),s(4:7,:),1)-1 ...                                     % Quaternion norm constraint
         -omega./k ...                                                      % Monotony of the time law
         R].';                                                              % Graspling inequality

    % Equality constraints
    ceq = [SF(4:6)+cross(params(7:9).', SF(1:3))-params(10:12).'; ...       % Velocity of the graspling fiture
          cos(tau(1,end))-cos(params(27)); ...                              % Multi-revolution anomaly constraint
          sin(tau(1,end))-sin(params(27)); ...                              % Multi-revolution anomaly constraint
          ];
end