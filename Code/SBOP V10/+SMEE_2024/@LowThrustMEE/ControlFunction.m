%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, S)
    % Constants
    mu = params(1);     % Gravitational parameter

    % Preallocation 
    u = zeros(3, size(S,2));
    old_u = u(3,:);

    % Auxiliary variables
    l = S(6,:);
    w = 1 + S(2,:) .* cos(l) + S(3,:) .* sin(l);
    delta = sqrt( S(1,:) / mu );
    dthetak = sqrt( mu * S(1,:) ) .* (w ./ S(1,:)).^2;
    s = 1 + S(4,:).^2 + S(5,:).^2;
    k = S(4,:) .* sin(l) - S(5,:) .* cos(l);
    dthetau = delta .* k ./ w;
    beta = delta .* s ./ (2 * w);
    
    % Linear terms of the equations of motion
    a = S(7:12,:);          % Inertial acceleration field

    % Normal component
    eps = 1E-10;            % Convergence tolerance
    iter = 1;               % Initial iteration
    max_iter = 100;         % Maximum number of iterations
    GoOn = true;            % Convergence flag
    
    % Iterations
    C = sqrt( a(4,:).^2 + a(5,:).^2 ) ./ beta;
    A = sign(a(4,:) ./ cos(l)); 
    dtheta = dthetak;
    u(3,:) = 0 * C .* A .* sign( dtheta );

    % Tangential component
    alpha = 2 * S(1,:) .* delta ./ w; 
    u(2,:) = a(1,:) ./ alpha;

    % Radial component
    for i = 1:size(S,2)
        B = OrbitalDynamics.MEE_matrix(mu, l(i), S(1:5,i));

        a_term = a(2,i) / delta(i) - B(2,2:3) * u(2:3,i); 
        b_term = a(3,i) / delta(i) - B(3,2:3) * u(2:3,i);
        u(1,i) = ( a_term )^2 + ( b_term )^2;

        Delta(1,1) = a_term / sin( l(i) );
        Delta(2,1) = b_term / cos( l(i) );
        
        u(1,i) = sqrt( u(1,i) );% * sign( Delta(1,1) ) * ( sign(Delta(1,1)) == sign(Delta(2,1)) );
    end
end