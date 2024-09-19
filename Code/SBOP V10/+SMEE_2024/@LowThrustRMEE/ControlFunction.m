%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, S)
    % Constants
    mu = params(1);     % Gravitational parameter

    % Preallocation 
    u = zeros(3, size(S,2));
%     old_u = u(3,:);

    % Auxiliary variables
    l = t(1,:);
    sin_l = sin(l);
    cos_l = cos(l);
    w = 1 + S(2,:) .* cos_l + S(3,:) .* sin_l;
    delta = sqrt( S(1,:) / mu );
    r_h = delta ./ w;
    dthetak = sqrt( mu * S(1,:) ) .* (w ./ S(1,:)).^2;
    s = 1 + S(4,:).^2 + S(5,:).^2;
    k = S(4,:) .* sin_l - S(5,:) .* cos_l;
    dthetau = r_h .* k;
    beta = r_h .* s / 2;
    
    % Linear terms of the equations of motion
    a = S(6:10,:);          % Inertial acceleration field

    % Normal component
%     idx_c = cos_l ~= 0;
%     Delta(2,idx_c) = a(4,idx_c) .* w ./ +cos_l(idx_c);
%     Delta(2,~idx_c) = zeros(1,sum(~idx_c));
% 
%     idx_s = sin_l ~= 0;
%     Delta(1,idx_s) = a(5,idx_s) .* w ./ +sin_l(idx_s);
%     Delta(1,~idx_s) = zeros(1,sum(~idx_s));
% 
%     eps = 1E-7;
%     sign_delta = tanh(Delta / eps);
%     A = sign_delta(1,:) + sign_delta(2,:);
%     A = A / 2 .* sign_delta(1,:);

    C = dot(a(4:5,:), a(4:5,:), 1) ./ beta.^2;
    
    ac = 1 - C .* dthetau.^2; 
    bc = - 2 * C .* dthetak .* dthetau; 
    cc = - C .* dthetak.^2;
    ah(1,:) = ( -bc + sqrt(bc.^2 - 4 * ac .* cc) ) ./ (2 * ac);
    ah(2,:) = ( -bc - sqrt(bc.^2 - 4 * ac .* cc) ) ./ (2 * ac);

    u(3,:) = min(ah, [], 1);
    dtheta = dthetak + dthetau .* u(3,:);

    % Tangential component
    alpha = 2 * S(1,:) .* r_h; 
    u(2,:) = a(1,:) .* dtheta ./ alpha;

    % Radial component
    ab_term = a(2:3,:) .* dtheta;

    for i = 1:size(S,2)
        B = OrbitalDynamics.MEE_matrix(mu, l(i), S(1:5,i));
        ab_term(:,i) = ab_term(:,i) - B(2:3,2:3) * u(2:3,i);
    end

%     Delta(1,idx_s) = ab_term(1,idx_s) ./ +sin_l(idx_s);
%     Delta(1,~idx_s) = zeros(1,sum(~idx_s));

%     Delta(2,idx_c) = ab_term(2,idx_c) ./ -cos_l(idx_c);
%     Delta(2,~idx_c) = zeros(1,sum(~idx_c));
% 
%     sign_delta = tanh(Delta / eps);
%     A = sign_delta(1,:) + sign_delta(2,:);
%     A = A / 2 .* sign_delta(1,:);

    u(1,:) = sqrt( dot(ab_term, ab_term, 1) ./ delta.^2 );
%     u(1,:) = u(1,:) .* A;
end