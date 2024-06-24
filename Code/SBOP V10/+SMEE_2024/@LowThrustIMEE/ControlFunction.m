%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, S)
    % Constants
    eps = 1E-5;
    mu = params(1);     % Gravitational parameter

    % Preallocation 
    u = zeros(3, size(S,2));

    % Auxiliary variables
    l = t(1,:) + S(6,:);
    w = 1 + S(2,:) .* cos(l) + S(3,:) .* sin(l);
    s = 1 + S(4,:).^2 + S(5,:).^2;
    delta = sqrt( S(1,:) / mu );
    k = S(4,:) .* sin(l) - S(5,:) .* cos(l);
    dtheta = sqrt(mu * S(1,:)) .* (w ./ S(1,:)).^2;
    
    % Linear terms of the equations of motion
    a = S(7:12,:);                     % Inertial acceleration field

    % Tangential component
    alpha = 2 * S(1,:) .* delta ./ (w .* dtheta); 
    u(2,:) = a(1,:) ./ alpha;

    % Normal component
    beta = ( k.^2 + s.^2 / 2 ) .* delta.^2 ./ (w .* dtheta).^2;
    u(3,:) = sqrt( a(4,:).^2 + a(5,:).^2 + a(6,:).^2 ./ beta );

    Delta(1,:) = a(4,:) ./ cos(l);
    Delta(2,:) = a(5,:) ./ sin(l);
    
    u(3,:) = u(3,:) .* sign( Delta(1,:) ) .* ( sign(Delta(1,:)) == sign(Delta(2,:)) );

    % Radial component
    for i = 1:size(S,2)
        B = OrbitalDynamics.MEE_matrix(mu, l(i), S(1:5,i));
        u(1,i) = ( a(2,i) .* dtheta(i) ./ delta(i) - B(2,2:3) * u(2:3,i))^2 + ( a(3,i) .* dtheta(i) ./ delta(i) - B(3,2:3) * u(2:3,i) )^2;

        Delta(1,1) = ( a(2,i) .* dtheta(i) ./ delta(i) - B(2,2:3) * u(2:3,i) ) / sin(l(i));
        Delta(2,1) = ( a(3,i) .* dtheta(i) ./ delta(i) - B(3,2:3) * u(2:3,i) ) / cos(l(i));
        u(1,i) = sqrt( u(1,i) ) * sign( Delta(1,1) ) * ( sign(Delta(1,1)) == sign(Delta(2,1)) );
    end
end