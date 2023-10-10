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

    % Compute the angular velocity
%     I = reshape(params(13:21), 3, 3);      % Inertia tensor of the chaser

    % Angular velocity
%     omega = 2 * [s(4,:).*s(11,:)+s(3,:).*s(12,:)-s(2,:).*s(13,:)-s(1,:).*s(14,:); ...
%                 -s(3,:).*s(11,:)+s(4,:).*s(12,:)+s(1,:).*s(13,:)-s(2,:).*s(14,:); ...
%                  s(2,:).*s(11,:)-s(1,:).*s(12,:)+s(4,:).*s(13,:)-s(3,:).*s(14,:)];
% 
%     % Angular acceleration 
%     alpha = 2 * [s(4,:).*s(18,:)+s(3,:).*s(19,:)-s(2,:).*s(20,:)-s(1,:).*s(21,:); ...
%                 -s(3,:).*s(18,:)+s(4,:).*s(19,:)+s(1,:).*s(20,:)-s(2,:).*s(21,:); ...
%                  s(2,:).*s(18,:)-s(1,:).*s(19,:)+s(4,:).*s(20,:)-s(3,:).*s(21,:)];
% 
%     Omega = [zeros(1,size(tau,2)); -Omega; zeros(1,size(tau,2))];
%     Alpha = [zeros(1,size(tau,2)); -Alpha; zeros(1,size(tau,2))];
%     for i = 1:size(tau,2)
%         aux = QuaternionAlgebra.RotateVector(s(4:7,i), Alpha(:,i));
%         Alpha(:,i) = aux(1:3,1);
%         aux = QuaternionAlgebra.RotateVector(s(4:7,i), Omega(:,i));
%         Omega(:,i) = aux(1:3,1);
%     end
%     omega_t = omega + Omega(1:3,:);
% 
%     for i = 1:length(tau)
%         u(4:6,i) = I * (alpha(:,i) + Alpha(:,i) - cross(omega(:,i), Omega(:,i))) + cross(omega_t(:,i), I*omega_t(:,i));           % Euler equations
%     end

    % Dimensionalization
    u(1:3,:) = u(1:3,:) .* Omega.^2 + Omega.*Alpha.*s(4:6,:);
%     u(4:6,:) = u(4:6,:) .* (rho.^2 .* k).^2 + (rho.^2 .* k).*(-2*k*params(24)*sin(tau(1,:).*rho)).*omega;
end