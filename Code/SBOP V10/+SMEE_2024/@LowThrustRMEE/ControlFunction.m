%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, S)
    % Constants
    mu = params(1);     % Gravitational parameter

    % Preallocation 
    u = zeros(3, size(S,2));

    % Auxiliary variables
    w = 1 + S(2,:) .* cos(t(1,:)) + S(3,:) .* sin(t(1,:));
    s = 1 + S(4,:).^2 + S(5,:).^2;
    delta = sqrt( S(1,:) / mu);
    
    % Linear terms of the equations of motion
    a = S(6:10,:);                     % Inertial acceleration field
    f = zeros(5,size(S,2));            % Acceleration vector

    % Keplerian term
    f(6,:) = sqrt(mu * S(1,:)) .* (w ./ S(1,:)).^2;

    % Tangential component
    alpha = 2 * S(1,:) ./ w .* delta; 
    u(2,:) = (a(1,:) - f(1,:)) ./ alpha;

    % Normal component
    beta = delta .* s ./ (2*w);
    u(3,:) = sqrt( a(4,:).^2 + a(5,:).^2 ) ./ beta;
    u(3,:) = u(3,:) .* sign( a(4,:) / cos(t(1,:)) );

    % Radial component
    for i = 1:size(S,2)
        B = OrbitalDynamics.MEE_matrix(mu, S(6,i), S(1:5,i));
        u(1,i) = ( a(2,i) - B(2,2:3) * u(2:3,i))^2 + ( a(3,i) - B(3,2:3) * u(2:3,i) )^2;
    end

    u(1,:) = sqrt( u(1,:) ) ./ delta .* sign( a(2,:) / sin(t(1,:)) );
end