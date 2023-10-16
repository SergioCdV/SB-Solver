%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Compute the auxiliary variables
    l = 2*s(1,:)-s(1,:).^2-s(2,:).^2;

    % Preallocation
    u = zeros(3,size(t,2));                            % Control vector
    u(2,:) = 0.5 * (s(9,:) +  s(2,:).*s(1,:).^3);      % Y acceleration
    u(1,:) = s(10,:)-(s(1,:)-1).*s(1,:).^3-u(2,:);     % X acceleration

    u(1,:) = u(1,:) ./ (s(3,:) .* l).^2;
    u(2,:) = u(2,:) ./ (s(3,:) .* l).^2;

    % Compute the complete osculating LVLH quaternion and angular velocity
    epsilon = s(4:7,:);

    omega = 2 * [epsilon(4,:).*s(12,:)+epsilon(3,:).*s(13,:)-epsilon(2,:).*s(14,:)-epsilon(1,:).*s(15,:); ...
                -epsilon(3,:).*s(12,:)+epsilon(4,:).*s(13,:)+epsilon(1,:).*s(14,:)-epsilon(2,:).*s(15,:); ...
                 epsilon(2,:).*s(12,:)-epsilon(1,:).*s(13,:)+epsilon(4,:).*s(14,:)-epsilon(3,:).*s(15,:)];

    % Normal acceleration
    u(3,:) = sqrt(dot(omega, omega, 1)) .* sign( s(12,:) ./ cos(s(8,:)) );

    % Dimensioning 
    u(3,:) = u(3,:) ./ (s(1,:) .* l).^2;
end