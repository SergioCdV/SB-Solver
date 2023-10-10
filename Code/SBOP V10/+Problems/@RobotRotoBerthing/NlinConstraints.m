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
%          dot(u(4:6,:),u(4:6,:),1)-params(22)^2 ...                       % Constraint on the torque magnitude (second order cone];                                                    
% 
%     omega_f = 2 * [s(4,end).*s(11,end)+s(3,end).*s(12,end)-s(2,end).*s(13,end)-s(1,end).*s(14,end); ...
%                 -s(3,end).*s(11,end)+s(4,end).*s(12,end)+s(1,end).*s(13,end)-s(2,end).*s(14,end); ...
%                  s(2,end).*s(11,end)-s(1,end).*s(12,end)+s(4,end).*s(13,end)-s(3,end).*s(14,end)];

    % Equality constraints
%     qc = QuaternionAlgebra.MPR2Quat(1, 1, params(7:9).', true);
%     dq = QuaternionAlgebra.right_isoclinic(s(4:7,end)) * QuaternionAlgebra.quaternion_inverse(qc);

    ceq = [cos(tau(1,end))-cos(params(7)); ...                               % Multi-revolution anomaly constraint
           sin(tau(1,end))-sin(params(7)); ...                               % Multi-revolution anomaly constraint
          % dq-[zeros(3,1);1]; ...                                           % Relative attitude constraint
          % omega_f-params(10:12).'; ...                                     % Relative angular velocity constraint
          % dot(s(4:7,:),s(4:7,:),1).'-1; ...                                % Quaternion norm constraint
          ];
end