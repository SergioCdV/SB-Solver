%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Orbit motion constants
%     rho = 1 + params(5) * cos(params(18:end));                     % Transformation parameter
%     k = params(4)^2/params(6)^3;                                   % Keplerian constant
%     Omega = rho.^2 .* k;                                           % True anomaly angular velocity [rad/s]
%     Alpha = -2 * k * params(5) * sin(params(18:end)) .* rho;       % True anomaly acceleration [rad/s^2]

    % Compute the angular velocity
    I = reshape(params(5:13), 3, 3);      % Inertia tensor of the chaser

    % Angular velocity
    omega = 2 * [s(4,:).*s(5,:)+s(3,:).*s(6,:)-s(2,:).*s(7,:)-s(1,:).*s(8,:); ...
                -s(3,:).*s(5,:)+s(4,:).*s(6,:)+s(1,:).*s(7,:)-s(2,:).*s(8,:); ...
                 s(2,:).*s(5,:)-s(1,:).*s(6,:)+s(4,:).*s(7,:)-s(3,:).*s(8,:)];

    % Angular acceleration 
    alpha = 2 * [s(4,:).*s(9,:)+s(3,:).*s(10,:)-s(2,:).*s(11,:)-s(1,:).*s(12,:); ...
                -s(3,:).*s(9,:)+s(4,:).*s(10,:)+s(1,:).*s(11,:)-s(2,:).*s(12,:); ...
                 s(2,:).*s(9,:)-s(1,:).*s(10,:)+s(4,:).*s(11,:)-s(3,:).*s(12,:)];

%     Omega = [zeros(1,size(tau,2)); -Omega; zeros(1,size(tau,2))];
%     Alpha = [zeros(1,size(tau,2)); -Alpha; zeros(1,size(tau,2))];
%     for i = 1:size(tau,2)
%         aux = QuaternionAlgebra.RotateVector(s(4:7,i), Alpha(:,i));
%         Alpha(:,i) = aux(1:3,1);
%         aux = QuaternionAlgebra.RotateVector(s(4:7,i), Omega(:,i));
%         Omega(:,i) = aux(1:3,1);
%     end
%     omega_t = omega + Omega(1:3,:);

    % Euler equations
    u = zeros(3,size(tau,2));
    for i = 1:length(tau)
        u(1:3,i) = I * alpha(:,i) + cross( omega(:,i), I*omega(:,i) );           
    end
end