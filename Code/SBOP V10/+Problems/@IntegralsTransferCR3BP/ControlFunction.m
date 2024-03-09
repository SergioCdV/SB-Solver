%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Constants 
    omega_p = params(3);         % In-plane frequency
    omega_v = params(4);         % Out-of-plane frequency
    k = params(5);               % In-plane constraint

    % Initial angles 
    psi = params(7);             % Out-of-plane angle
    phi = params(8);             % In-plane angle

    t = tau(1,:);
    
    % Compute the dipoles
    f = [-2 * (omega_p - k) * sin(omega_p * t + phi) .* s(3,:); ...
         +2 * omega_v .* cos(omega_v * t + psi) .* s(4,:)];

    % Compute the control vector as a dynamics residual
    u([1 3],:) = [s(5,:) .* cos(omega_p * t + phi); s(6,:) .* sin(omega_v * t + psi)] + f;
    u(1,:) = -u(1,:);

    % Final constraint 
    u(2,:) = k * s(5,:) .* sin(omega_p * t + phi) + 2 * (omega_p * k - 1) * cos(omega_p * t + phi) .* s(3,:); 
end