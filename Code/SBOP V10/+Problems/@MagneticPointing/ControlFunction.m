%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Parameters
    COE = params(2:7);                      % Orbit elements
    I = reshape(params(8:16), [3 3]);       % Inertia matrix

    % Orbital dynamics 
    v = tau(1,:);                                           % True anomaly space
    h = sqrt(params(1) * COE(1) * (1-COE(2)^2));            % Angular momentum 
    r = COE(1) * (1-COE(2)^2) ./ (1 + COE(2)*cos(v));       % Spacecraft position vector
    dtheta = h./r.^2;                                       % Angular velocity of the anomaly

    % Angular velocity
    omega = 2 * [s(4,:).*s(5,:)+s(3,:).*s(6,:)-s(2,:).*s(7,:)-s(1,:).*s(8,:); ...
                -s(3,:).*s(5,:)+s(4,:).*s(6,:)+s(1,:).*s(7,:)-s(2,:).*s(8,:); ...
                 s(2,:).*s(5,:)-s(1,:).*s(6,:)+s(4,:).*s(7,:)-s(3,:).*s(8,:)];

    % Angular acceleration 
    alpha = 2 * [s(4,:).*s(9,:)+s(3,:).*s(10,:)-s(2,:).*s(11,:)-s(1,:).*s(12,:); ...
                -s(3,:).*s(9,:)+s(4,:).*s(10,:)+s(1,:).*s(11,:)-s(2,:).*s(12,:); ...
                 s(2,:).*s(9,:)-s(1,:).*s(10,:)+s(4,:).*s(11,:)-s(3,:).*s(12,:)];

    % Compute the control vector as a dynamics residual
    u = zeros(3,size(s,2));
    for i = 1:size(s,2)
        u(:,i) = (dtheta(i) * I * alpha(:,i) + cross( omega(:,i), I * omega(:,i) )) / beta(1);
    end

%     u = beta(end) * u;
end