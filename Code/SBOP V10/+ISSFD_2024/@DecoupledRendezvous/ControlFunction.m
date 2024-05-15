%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Constants
    k = params(2)^2 / params(4)^3;                % True anomaly angular velocity
    rho = 1 + params(3) * cos(tau(1,:));          % Transformation parameter
    drho = -params(3) .* sin(tau(1,:));           % Derivative of the transformation parameter
    olvlh = k * rho.^2;                           % Angular velocity of the target's LVLH frame

    % Compute the control vector as a dynamics residual (linear acceleration, TH relative motion model)
    u(1:3,:) = s(13:15,:) - [2 * s(9,:); -s(2,:); 3 * s(3,:) ./ rho - 2 * s(7,:)];

    % Dimensionalisation 
    u(1:3,:) = u(1:3,:) .* olvlh.^2 ./ rho;

    % Attitude control 
%     I = reshape(params(9:17), [3 3]);       % Inertia tensor of the chaser
% 
%     % Shadow transformation
%     sigma = s(4:6,:);
%     idx = dot(sigma, sigma, 1) > 1;
%     sigma(:,idx) = -sigma(:,idx) ./ dot(sigma(:,idx), sigma(:,idx), 1);
% 
%     % Pre-allocation
%     q = [sigma; -ones(1,size(tau,2))];      % Modified MRPs
%     q_squared = dot( q, q ).^2;             % Dot product of the associated quaternions
%     alpha = zeros(3,size(tau,2));           % Angular velocity of the target
%         
%     % Angular velocity
%     dsigma =  s(10:12,:) .* olvlh;
%     ddsigma = s(16:18,:) .* olvlh.^2 + (- 2 * k .* rho .* drho) .* dsigma;
% 
%     omega = dsigma ./ q_squared;
%     
%     for i = 1:size(tau,2)
%         B = QuaternionAlgebra.Quat2Matrix( q(:,i) ).';
%         omega(:,i) = 4 * B * omega(:,i);
% 
%         dB = QuaternionAlgebra.dBmatrix( sigma(:,i), omega(:,i) );
%         alpha(:,i) = 4 * B * ( ddsigma(:,i) - 0.25 * dB * omega(:,i) );
%     end
% 
%     alpha = alpha ./ q_squared;                                                 % Angular acceleration
%     alpha = (I * alpha) + 2 * 0 * olvlh .* rho .* drho .* (I * omega);              % Derivative of the angular momentum
%     u(4:6,:) = alpha + cross( omega, I * omega );                               % Torque on the spacecraft 
    u(4:6,:) = zeros(3,size(tau,2));
end