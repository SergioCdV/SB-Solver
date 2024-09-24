%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Preallocation
    u = zeros(3,size(t,2));                        % Control vector

    % Compute the auxiliary variables
    sigma = t(1,:); 
    cos_sigma = cos(sigma); 
    sin_sigma = sin(sigma);
    l = 1 + s(1,:) .* cos_sigma + s(2,:) .* sin_sigma;
    gamma = l.^3 .* s(3,:).^4;

    % Accelerations 
    a = s(8:14,:);

    % Normal acceleration
    epsilon = s(4:7,:);         % Attitude quaternion

    omega = 2 * [epsilon(4,:) .* a(4,:) + epsilon(3,:) .* a(5,:) - epsilon(2,:) .* a(6,:) - epsilon(1,:) .* a(7,:); ...
                 epsilon(4,:) .* a(5,:) + epsilon(1,:) .* a(6,:) - epsilon(3,:) .* a(4,:) - epsilon(2,:) .* a(7,:); ...
                 epsilon(4,:) .* a(6,:) - epsilon(2,:) .* a(4,:) - epsilon(1,:) .* a(5,:) - epsilon(3,:) .* a(7,:)];
    
    omega = omega ./ sqrt( dot(epsilon, epsilon, 1) );

    u(3,:) = sqrt( dot(omega, omega, 1) ) .* gamma;

%     Delta = zeros(1, size(u,2));
%     idx_c = cos_sigma ~= 0; 
%     idx_s = sin_sigma ~= 0;
%     idx = idx_c & idx_s; 
%     Delta(idx) = sign( omega(1,idx) ./ cos_sigma(idx) );
%     % .* Delta;

    % Tangential acceleration
    u(2,:) = -a(3,:) .* (s(3,:).^3 .* l.^3);      

    % Radial acceleration
    B = [s(1,:) + (1 + l) .* cos_sigma; s(2,:) + (1 + l) .* sin_sigma];
    ab_term = a(1:2,:) - B .* u(2,:) ./ gamma ;
    u(1,:) = sqrt( dot(ab_term, ab_term, 1) ) .* gamma ./ l;

%     Delta = zeros(1, size(u,2));
%     idx_c = cos_sigma ~= 0; 
%     idx_s = sin_sigma ~= 0;
%     idx = idx_c & idx_s;
%     Delta(idx) = sign( a(1,idx) ./ sin_sigma(idx) );
end