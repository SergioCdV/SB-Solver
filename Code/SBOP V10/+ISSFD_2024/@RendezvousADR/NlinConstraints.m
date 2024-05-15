%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % TH transformation matri
    k = params(2)^2 / params(4)^3;                                      % True anomaly angular velocity
    rho = 1 + params(3) * cos(tau(1,:));                                % Transformation parameter
    drho =  - params(3) * sin(tau(1,:));                                % Derivative of the transformation
    omega = k * rho.^2;                                                 % Angular velocity of the target's LVLH frame

    % Angular velocity of the chaser
    sigma = s(4:6,:);
    idx = dot(sigma, sigma, 1) > 1;
    sigma(:,idx) = -sigma(:,idx) ./ dot(sigma(:,idx), sigma(:,idx), 1);

    q = [sigma; -ones(1,size(tau,2))];                    % Modified MRPs
    Omega = s(10:12,:) .* omega ./ dot(q,q,1);            % Angular velocity of the chaser

    % Linear velocity of the chaser
    v = zeros(3,size(tau,2));

    for i = 1:size(tau,2)
        A = [rho(i) * eye(3) zeros(3); drho(i) * eye(3) eye(3)/(rho(i) * omega(i))];    % TH transformation matrix
        v_aux = A \ [s(1:3,i); s(7:9,i)];
        v(:,i) = v_aux(4:6,1);

        B = QuaternionAlgebra.Quat2Matrix( q(:,i) );
        Omega(:,i) = 4 * B.' * Omega(:,i);
    end
    
    % Final position in physical space    
    rb = s(1:3,:) ./ rho;
    rf = rb(1:3,end); 

    % LOS constraint
    cx = params(27); 
    cz = params(28); 
    xmin = params(29); 
    zmin = params(30);

    Alos = [0 -1 0; cx -1 0; -cx -1 0; 0 -1 cz; 0 -1 -cz];
    Clos = [0; cx * xmin; cx * xmin; cz * zmin; cz * zmin]; 

    qc = QuaternionAlgebra.MPR2Quat(1, 1, sigma, true);                     % Quaternion between the two body frames 
    inv_qc = QuaternionAlgebra.quaternion_inverse(qc);                      % Rotation from chaser body frame to LVLH axes

    qt = reshape(params(36:35 + 4 * (params(35)+1)), 4, params(35)+1);

    for i = size(tau,2):size(tau,2)
        q = QuaternionAlgebra.right_isoclinic(qt(:,i)) * inv_qc(:,i);       % Quaternion from chaser body frame to target body frame
        rb(:,i) = QuaternionAlgebra.RotateVector(q, rb(:,i));               % Rotation to target axes
    end

    LOS = Alos * rb - Clos;

    % LOS pointing constraint
    dock_chaser = reshape(params(18:20), 3, 1);                     % Docking port of the chaser in the chaser's body frame
    dock_chaser = QuaternionAlgebra.RotateVector(q, dock_chaser);   % Docking port of the chaser in the target's body frame
    target_chaser = reshape(params(21:23), 3, 1);                   % Docking port of the target in the target's body frame

    LOS_angle = params(31) - dot(dock_chaser, -target_chaser); 

    % Angular velocity residual 
    omega_f = QuaternionAlgebra.RotateVector(q, Omega(:,end));   % Final chaser's body frame angular velocity with respect to LVLH in the target's body frame
    omega_res = omega_f - params(32:34).';                       % Residual to the target's angular velocity

    % Inequalities
    c = [
            dot(u(1:3,:), u(1:3,:), 1) - params(7)^2 ...         % Constraint on the force magnitude (second order cone)
            reshape(+v - params(24), 1, []) ...                  % Maximum linear velocity of the chaser
            reshape(-v - params(24), 1, []) ...                  % Maximum linear velocity of the chaser
            reshape(+u(4:6,:) - params(8), 1, []) ...            % Constraint on the torque magnitude (infinty norm)
            reshape(-u(4:6,:) - params(8), 1, [])...             % Constraint on the torque magnitude (infinty norm)
            reshape(+Omega - params(25), 1, []) ...              % Maximum angular velocity of the chaser
            reshape(-Omega - params(25), 1, []) ...              % Maximum angular velocity of the chaser
            LOS_angle ...                                        % LOS angle constraint
            dot(rf, rf, 1) - params(26)^2 ...                    % Graspling constraint
            reshape(LOS, 1, []) ...                              % LOS constraint
%             reshape(+omega_res - 1E-4, 1, []) ...                % Angular velocity sync constraint
%             reshape(-omega_res - 1E-4, 1, []) ...                % Angular velocity sync constraint
        ];      

    % Equality constraints
    ceq = [];
end