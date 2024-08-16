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
    sin_l = sin(l);
    cos_l = cos(l);
    w = 1 + S(2,:) .* cos_l + S(3,:) .* sin_l;
    delta = sqrt( S(1,:) / mu );
    r_h = delta ./ w;
    dtheta = sqrt(mu * S(1,:)) .* (w ./ S(1,:)).^2;
    k = S(5,:) .* cos_l - S(4,:) .* sin_l;
    
    % Linear terms of the equations of motion
    a = S(7:12,:);                       % Inertial acceleration field

    % Regularization 
    a = a .* dtheta;                    

    % Tangential component
    alpha = 2 * S(1,:) .* r_h; 
    u(2,:) = a(1,:) ./ alpha;

    % Normal componen.t
    beta = ( r_h.^2 ./ (16 * (1 - sigma_norm).^2) ) .* ( (1 - sigma_norm.^2).^2 + k.^2 .* (64 + 4 * (1 + sigma_norm).^2) );
    u(3,:) = sqrt( dot( a(4:6,:), a(4:6,:), 1 ) ./ beta );

    A = (1 - S(4,:).^2 + S(5,:).^2) .* cos_l - 2 * S(4,:) .* S(5,:) .* sin_l;
    idx_c = A ~= 0;
    Delta(2,idx_c) = a(4,idx_c) ./ +A(idx_c);
    Delta(2,~idx_c) = zeros(1,sum(~idx_c));
    
    B = (1 + S(4,:).^2 - S(5,:).^2) .* sin_l - 2 * S(4,:) .* S(5,:) .* cos_l;
    idx_s = B ~= 0;
    Delta(1,idx_s) = a(5,idx_s) ./ +B(idx_s);
    Delta(1,~idx_s) = zeros(1,sum(~idx_s));

    A = sign( Delta(1,:) ) .* ( sign(Delta(1,:)) == sign(Delta(2,:)) );
    u(3,:) = u(3,:) .* A;

    % Radial component
    ab_term = a(2:3,:);

    for i = 1:size(S,2)
        B = OrbitalDynamics.SMEE_matrix(mu, l(i), S(1:5,i));
        ab_term(:,i) = ab_term(:,i) - B(2:3,2:3) * u(2:3,i); 
    end

    ab_term = ab_term ./ delta;

    idx_s = sin_l ~= 0;
    Delta(1,idx_s) = ab_term(1,idx_s) ./ +sin_l(idx_s);
    Delta(1,~idx_s) = zeros(1,sum(~idx_s));

    idx_c = cos_l ~= 0;
    Delta(2,idx_c) = ab_term(2,idx_c) ./ -cos_l(idx_c);
    Delta(2,~idx_c) = zeros(1,sum(~idx_c));

    u(1,:) = sqrt( dot(ab_term, ab_term, 1) );
    u(1,:) = u(1,:) .* sign( Delta(1,:) ) .* ( sign(Delta(1,:)) == sign(Delta(2,:)) );
end