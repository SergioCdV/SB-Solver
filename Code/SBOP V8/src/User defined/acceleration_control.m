%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the 9xm state vector 
%         - scalar tf, the final time of flight

% Outputs: - vector u, the nondimensional 3xm control vector
%          - vector dv, the inertial velocity field 
%          - vector f, the dynamics vector field

function [u, dv, f] = acceleration_control(mu, C, tf)

    % Set up value iteration 
    

    % Linear terms of the equations of motion
    c = tf;                                                                         % Normalizing factor
    a = C(9:12,:);                                                                  % Inertial acceleration field
    f = zeros(size(a));      % Acceleration vector
    dv = C(5:8,:);                                                                  % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = a-f;
    for i = 1:size(C,2)
        L = KS_matrix(C(1:4,i));
        u(:,i) = L*u(:,i);
    end
end
