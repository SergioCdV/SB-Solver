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
    % Radial distance 
    r = dot(C(1:4,:),C(1:4,:),1); 

    % Energy evolution 
    E = -4*(tf^2*mu/2-dot(C(5:8,:),C(5:8,:),1))./r;

    % Linear terms of the equations of motion
    a = C(9:12,:);                                   % Inertial acceleration field
    f = -(E./2).*C(1:4,:);                           % Acceleration vector
    dv = C(5:8,:);                                   % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = 2*(a-f)./r;
end
