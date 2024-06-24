%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, S)
    % Constants
    mu = params(1);     % Gravitational parameter

    % Preallocation 
    u = zeros(3, size(S,2));

    % Shadow set transformation 
    sigma_norm = dot(S(4:5,:), S(4:5,:), 1);
    test = sigma_norm > 1;                             % Shadow set switching function
    S(4:5,test) = - S(4:5,test) ./ sigma_norm(test);   % Shadow set 
    sigma_norm = dot(S(4:5,:), S(4:5,:), 1);           % New norm

    % Auxiliary variables
    l = t(1,:) + S(6,:);
    w = 1 + S(2,:) .* cos(l) + S(3,:) .* sin(l);
    delta = sqrt( S(1,:) / mu );
    k = S(5,:) .* cos(l) - S(4,:) .* sin(l);
    dtheta = sqrt(mu * S(1,:)) .* (w ./ S(1,:)).^2;
    
    % Linear terms of the equations of motion
    a = S(7:12,:);                     % Inertial acceleration field

    % Tangential component
    alpha = 2 * S(1,:) .* delta ./ (w .* dtheta); 
    u(2,:) = a(1,:) ./ alpha;

    % Normal component
    beta = (1 - sigma_norm).^2 + (64 * (1 + sigma_norm).^2 - 4) .* k.^2;
    u(3,:) = 4 * (1-sigma_norm.^2) .* delta .* dtheta .* sqrt( a(4,:).^2 + a(5,:).^2 + a(6,:).^2 ./ beta );
    u(3,:) = u(3,:) .* sign( a(6,:) .* dtheta / (delta .* (-k) ./ w) );

    % Radial component
    for i = 1:size(S,2)
        B = OrbitalDynamics.SMEE_matrix(mu, l(i), S(1:5,i));
        u(1,i) = ( a(2,i) .* dtheta(i) ./ delta(i) - B(2,2:3) * u(2:3,i))^2 + ( a(3,i) .* dtheta(i) ./ delta(i) - B(3,2:3) * u(2:3,i) )^2;
        u(1,i) = sqrt( u(1,i) ) .* sign( ( a(2,i) .* dtheta(i) ./ delta(i) - B(2,2:3) * u(2:3,i) ) / sin(l(i)) );
    end
end