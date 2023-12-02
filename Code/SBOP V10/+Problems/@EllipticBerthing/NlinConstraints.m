%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % TH transformation matrix
    k = params(4)^2 / params(6)^3;                                           % True anomaly angular velocity
    rho = 1 + params(5) * cos(tau(1,:));                                     % Transformation parameter
    kp =  - params(5) * sin(tau(1,end));                                     % Derivative of the transformation
    L = [rho(end) * eye(3) zeros(3); kp * eye(3) eye(3)/(k * rho(end))];     % TH transformation matrix
    SF = L \ s([1:3 7:9], end);       

    % Variables 
    Rho = s(1:3,:);             % Relative position
    sigma = s(4:6,:);           % MRP
    drho = s(7:9,:);            % Relative velocity
    dsigma = s(10:12,:);        % Derivative of the MRP

    % Integration of the target's attitude 
    q_t = params(37:40).';      % Final target attitude
    omega_t = params(41:43).';  % Final target angular velocity

    % Inequality constraints
    R1 = +dot(SF(1:3,1), SF(1:3,1)) - params(9)^2; 
    L =  -dot(SF(1:3,1), SF(1:3,1)) + params(10)^2;

    omega = zeros(3,length(tau));

    for i = 1:length(tau)
        q = [sigma(:,i); -1];
        B = QuaternionAlgebra.Quat2Matrix(q);
        omega(:,i) = 4 * B.' * dsigma(:,i) / dot(q,q)^2;          
    end

    q_c = QuaternionAlgebra.MPR2Quat(1, 1, sigma(:,end), true);
    q_e = QuaternionAlgebra.quaternion_inverse( q_c );
    r = params(31:33).';
    r = QuaternionAlgebra.RotateVector(q_e, r);

    q_inv = QuaternionAlgebra.quaternion_inverse( q_t(:,end) );
    b = params(34:36).';
    b = QuaternionAlgebra.RotateVector(q_inv, b);

    omega_r = omega(:,end) - QuaternionAlgebra.RotateVector( QuaternionAlgebra.right_isoclinic(q_c) * q_inv, omega_t(:,end) );

    theta_c = params(44);
    theta_e = params(45);

    z = Rho(:,end) / norm(Rho(:,end));

    c = [ 
          sqrt(dot(u(1:3,:), u(1:3,:), 1)) - params(3) ...          % Constraint on the force magnitude (second order cone)
          R1 ...                                                    % Keep-out sphere constraint
          L ...                                                     % Graspling reach
          dot(sigma, sigma)-1 ...                                   % MRP kinematic constraint
          +reshape(omega, 1, []) - params(21) ...                   % Infinity norm of the angular velocity
          -reshape(omega, 1, []) - params(21) ...                   % Infinity norm of the angular velocity
          +reshape(u(4:6,:), 1, []) - params(20) ...                % Attitude control authority
          -reshape(u(4:6,:), 1, []) - params(20) ...                % Attitude control authority
          cos(theta_c) - dot(z, r) ...                              % Corridor constraint
          cos(theta_e) + dot(r, b) ...                              % Cone constraint
        ];
    
    % Equality constraints
    ceq = [
            SF(4:6,end); ...                                        % Final relative veleocity
            omega_r                                                 % Final angular velocity                                                                
          ];
end