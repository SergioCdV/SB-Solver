%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    k = params(4)^2 / params(6)^3;                               % True anomaly angular velocity
    rho = 1 + params(5) * cos(tau(1,:));                         % Transformation parameter
    Omega = k .* rho.^2;                                         % Orbital motion [rad/s]
    dOmega = -2 * k .* rho .* (params(5) * sin(tau(1,:)));       % Derivative of the orbital motion

    % Compute the control vector as a dynamics residual (linear acceleration, TH relative motion model)
    u(1:3,:) = s(13:15,:) - [2 * s(9,:); -s(8,:); 3*s(3,:)./rho-2*s(7,:)];

    % Attitude control
    I = reshape(params(11:19), 3, 3);      % Inertia tensor of the chaser
    omega = zeros(3,length(tau));
    alpha = omega;

    for i = 1:length(tau)
        q = [s(4:6,i); -1];
        B = QuaternionAlgebra.Quat2Matrix(q);
        omega(:,i) = 4 * B.' * Omega(i) * s(10:12,i) / dot(q,q)^2;
        dB = dBmatrix(s(4:6,i), s(10:12,i));
        alpha(:,i) = 4 * B.' * (Omega(i) * s(16:18,i) + dOmega(i) * s(10:12,i) - 0.25 * dB * omega(:,i)) / dot(q,q)^2;
        u(4:6,i) = I * alpha(:,i) + cross( omega(:,i), I*omega(:,i) );           
    end

    % Dimensionalisation 
    u(1:3,:) = u(1:3,:) .* (k^2 * rho.^3);
    % u(4:6,i) = 
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