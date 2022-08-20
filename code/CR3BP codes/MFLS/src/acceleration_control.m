%% CR3BP Library %% 
% Sergio Cuevas del Valle
% Date: 18/08/22
% File: mfls_optimization.m 
% Issue: 0 
% Validated: 18/08/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - array C, the 9xm state vector
%         - scalar tf, the final time of flight

% Outputs: - vector u, the nondimensional 3xm control vector

function [u] = acceleration_control(wp, wv, kap, phi, psi, C, tf, tau) 
    % Compute the control vector as a residual
    C(3:4,:) = C(3:4,:)/tf; 
    C(5:6,:) = C(5:6,:)/tf^2;
    u(1,:) = -C(5,:).*cos(tf*tau*wp+phi)+2*wp*C(3,:).*sin(tf*tau*wp+phi);
    u(2,:) = kap*(C(5,:).*sin(tf*tau*wp+phi)+2*wp*C(3,:).*cos(tf*tau*wp+phi));
    u(3,:) = C(6,:).*sin(tf*tau*wv+psi)+2*wv*C(4,:).*cos(tf*tau*wv+psi);
end
