%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Pre-allocation
    sigma = s(1:3,:);
    dsigma = s(7:9,:);

    % Shadow transformation
    idx = dot(sigma, sigma, 1) > 1;
    sigma(:,idx) = -sigma(:,idx) ./ dot(sigma(:,idx), sigma(:,idx), 1);

    q = [sigma; -ones(1,size(tau,2))];      % Modified MRPs
    q_squared = dot( q, q ).^2;             % Dot product of the associated quaternions
  
    omega = dsigma ./ q_squared;
    
    for i = 1:size(tau,2)
        B = QuaternionAlgebra.Quat2Matrix( q(:,i) ).';
        omega(:,i) = 4 * B * omega(:,i);
    end

    % Inequality constraints
    res = u + cross(omega, s(4:6,:)) + s(10:12,:);     % Constraint between the derivative of the angular momentum and the control vector
    
    c = [
            dot(u,u,1) - params(3)^2 ...                    % Euclidean norm on the torque
            reshape(+s(4:6,:) - params(5), 1, []) ...       % Infinity norm on the RW
            reshape(-s(4:6,:) - params(5), 1, []) ...       % Infinity norm on the RW
            reshape(+omega - params(4), 1, []) ...          % Infinity norm on the angular velocity
            reshape(-omega - params(4), 1, []) ...          % Infinity norm on the angular velocity
            reshape(+res - 5E-7, 1, []) ...          % Infinity norm on the angular velocity
            reshape(-res - 5E-7, 1, []) ...          % Infinity norm on the angular velocity
        ]; 

    % Equality constraints
    ceq = [
%             res
          ];
end