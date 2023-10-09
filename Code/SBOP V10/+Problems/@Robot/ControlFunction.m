%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Compute the control vector as a dynamics residual (linear acceleration, TH relative motion model)
    Omega = params(23)^2 / params(25)^3;                                    % True orbit mean motion
    rho = 1 + params(24) * cos(tau(1,:));                                   % Transformation parameter
    u(1:3,:) = s(15:17,:) - [2 * s(3,:); -s(2,:); 2*s(3,:)./rho-2*s(1,:)];

    % Compute the angular velocity
    I = reshape(params(13:21), 3, 3);      % Inertia tensor of the chaser

    % Angular velocity
    omega = 2 * [s(4,:).*s(11,:)+s(3,:).*s(12,:)-s(2,:).*s(13,:)-s(1,:).*s(14,:); ...
                -s(3,:).*s(11,:)+s(4,:).*s(12,:)+s(1,:).*s(13,:)-s(2,:).*s(14,:); ...
                 s(2,:).*s(11,:)-s(1,:).*s(12,:)+s(4,:).*s(13,:)-s(3,:).*s(14,:)];

    % Angular acceleration 
    alpha = 2 * [s(4,:).*s(18,:)+s(3,:).*s(19,:)-s(2,:).*s(20,:)-s(1,:).*s(21,:); ...
                -s(3,:).*s(18,:)+s(4,:).*s(19,:)+s(1,:).*s(20,:)-s(2,:).*s(21,:); ...
                 s(2,:).*s(18,:)-s(1,:).*s(19,:)+s(4,:).*s(20,:)-s(3,:).*s(21,:)];

    for i = 1:length(tau)
        u(4:6,i) = I * alpha(:,i) + cross(omega(:,i), I*omega(:,i));           % Euler equations
    end

    % Dimensionalization
    u = u .* (Omega.*rho).^2;
end