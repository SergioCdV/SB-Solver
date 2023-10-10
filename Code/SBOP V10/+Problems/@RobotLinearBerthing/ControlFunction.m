%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    rho = 1 + params(5) * cos(tau(1,:));                             % Transformation parameter
    k = params(4)^2/params(6)^3;                                     % Keplerian constant
    Omega = rho.^2 .* k;                                             % True anomaly angular velocity [rad/s]
    Alpha = -2 * k * params(5) * sin(tau(1,:)) .* rho;               % True anomaly acceleration [rad/s^2]

    % Compute the control vector as a dynamics residual (linear acceleration, TH relative motion model)
    u(1:3,:) = s(7:9,:) - [2 * s(3,:); -s(2,:); 2*s(3,:)./rho-2*s(1,:)];

    % Dimensionalization
    u(1:3,:) = u(1:3,:) .* Omega.^2 + Omega.*Alpha.*s(4:6,:);
end