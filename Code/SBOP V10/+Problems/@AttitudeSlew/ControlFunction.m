%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Constants
    I = reshape(params(6:end), [3 3]);      % Inertia matrix

    % Pre-allocation
    omega = zeros(3, size(tau,2));          % Angular velocity
    alpha = omega;                          % Angular acceleration
    sigma = s(1:3,:);
    dsigma = s(7:9,:);
    ddsigma = s(13:15,:);

    % Shadow transformation
    idx = dot(sigma, sigma, 1) > 1;
    sigma(:,idx) = -sigma(:,idx) ./ dot(sigma(:,idx), sigma(:,idx), 1);

    q = [sigma; -ones(1,size(tau,2))];      % Modified MRPs
    q_squared = dot( q, q ).^2;             % Dot product of the associated quaternions
  
    omega = dsigma ./ q_squared;
    
    for i = 1:size(tau,2)
        B = QuaternionAlgebra.Quat2Matrix( q(:,i) ).';
        omega(:,i) = 4 * B * omega(:,i);

        dB = QuaternionAlgebra.dBmatrix( sigma(:,i), omega(:,i) );
        alpha(:,i) = 4 * B * ( ddsigma(:,i) - 0.25 * dB * omega(:,i) );
    end

    alpha = alpha ./ q_squared;                                                 % Angular acceleration
    
    % Torque solution
    u = I * alpha + cross(omega, I * omega);
end