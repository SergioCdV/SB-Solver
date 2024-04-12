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

    % Linear velocity of the chaser
    v = zeros(3,size(tau,2));
    for i = 1:size(tau,2)
        A = [rho(i) * eye(3) zeros(3); drho(i) * eye(3) eye(3)/(rho(i) * omega(i))];    % TH transformation matrix
        v_aux = A \ [s(1:3,i); s(7:9,i)];
        v(:,i) = v_aux(4:6,1);
    end

    % Angular velocity of the chaser
    sigma = s(4:6,:);
    idx = sqrt(dot(sigma, sigma, 1)) > 1;
    sigma(:,idx) = -sigma(:,idx) ./ dot(sigma(:,idx), sigma(:,idx), 1);

    q = [sigma; -ones(1,size(tau,2))];                    % Modified MRPs
    omega = s(10:12,:) .* omega ./ dot(q,q,1);            % Angular velocity of the chaser

    for i = 1:size(tau,2)
        B = QuaternionAlgebra.Quat2Matrix( q(:,i) );
        omega(:,i) = 4 * B.' * omega(:,i);
    end

    c = [
            dot(u(1:3,:), u(1:3,:), 1) - params(7)^2 ...         % Constraint on the force magnitude (second order cone)
            dot(u(4:6,:), u(4:6,:), 1) - params(8)^2 ...         % Constraint on the torque magnitude (second order cone)
            reshape(+v - params(24), 1, []) ...                  % Maximum linear velocity of the chaser
            reshape(-v - params(24), 1, []) ...                  % Maximum linear velocity of the chaser
%             reshape(+omega - params(25), 1, []) ...              % Maximum angular velocity of the chaser
%             reshape(-omega - params(25), 1, []) ...              % Maximum angular velocity of the chaser
        ];         

    % Equality constraints
    ceq = [];
end