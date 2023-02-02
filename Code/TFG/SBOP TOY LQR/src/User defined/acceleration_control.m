%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - array C, the 3xm state vector 
%         - scalar tf, the final time of flight
%         - vector t, the independent variable evolution

% Outputs: - vector u, the nondimensional 3xm control vector
%          - vector dv, the inertial velocity field 
%          - vector f, the dynamics vector field

function [u, dv, f] = acceleration_control(C, tf, t)
    % Compute the dynamics residual
    dv = C(2,:);           % Velocity field
    f = C(3,:);            % Newtonian dynamics
    u = C(3,:);            % Control acceleration
end
