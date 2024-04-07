%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Constants
    k = params(2)^2 / params(4)^3;                % True anomaly angular velocity
    rho = 1 + params(3) * cos(tau(1,:));          % Transformation parameter
    olvlh = k * rho.^2;                           % Angular velocity of the target's LVLH frame

    % Compute the control vector as a dynamics residual (linear acceleration, TH relative motion model)
    u(1:3,:) = s(13:15,:) - [2 * s(9,:); -s(2,:); 3 * s(3,:) ./ rho - 2 * s(7,:)];

    % Dimensionalisation 
    u(1:3,:) = u(1:3,:) .* olvlh.^2 ./ rho;

    % Attitude control 
%     I = reshape(params(9:17), [3 3]);       % Inertia tensor of the chaser

    % Shadow transformation
%     sigma = s(4:6,:);
%     idx = dot(sigma, sigma, 1) > 1;
%     sigma(:,idx) = -sigma(:,idx) ./ dot(sigma(:,idx), sigma(:,idx), 1);
% 
%     q = [sigma; -ones(1,size(tau,2))];      % Modified MRPs
%     omega = zeros(3,size(tau,2));           % Angular velocity of the chaser
%     alpha = omega;                          % Angular velocity of the target
% 
%     dsigma =  s(10:12,:) .* olvlh;
%     ddsigma = s(16:18,:) .* olvlh.^2 + 2 * k .* rho .* (-params(3) * sin(tau(1,:))) .* olvlh .* dsigma;
% 
%     for i = 1:size(tau,2)
%         B = QuaternionAlgebra.Quat2Matrix(q);
%         omega(:,i) = 4 * B.' * dsigma(:,i) / dot( q(:,i), q(:,i) )^2;
%         dB = dBmatrix(sigma(:,i), dsigma(:,i));
%         alpha(:,i) = 4 * B.' * ( ddsigma(:,i) - 0.25 * dB * omega(:,i) );
%     end
%     
%     u(4:6,:) = I * alpha ./ dot(q, q, 1).^2 + cross( omega, I * omega );  % Torque on the spacecraft 
end

%% Auxiliary function 
function [Bdot] = dBmatrix(q, dq)
    Bdot = zeros(3,3);
    s = -2 * dot(q, dq);
    Bdot(1,1) = s + 4 * (q(1) * dq(1));
    Bdot(1,2) = 2 * (-dq(3) + q(1) * dq(2) + dq(1) * q(2));
    Bdot(1,3) = 2 * ( dq(2) + q(1) * dq(3) + dq(1) * q(3));
    Bdot(2,1) = 2 * ( dq(3) + q(1) * dq(2) + dq(1) * q(2));
    Bdot(2,2) = s + 4 * (q(2) * dq(2));
    Bdot(2,3) = 2 * (-dq(1) + q(2) * dq(3) + dq(2) * q(3));
    Bdot(3,1) = 2 * (-dq(2) + q(1) * dq(3) + dq(1) * q(3));
    Bdot(3,2) = 2 * ( dq(1) + q(2) * dq(3) + dq(2) * q(3));
    Bdot(3,3) = s + 4 * (q(3) * dq(3));
end