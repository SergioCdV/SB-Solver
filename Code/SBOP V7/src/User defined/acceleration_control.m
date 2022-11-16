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
    % Compute the radius vector
    r = sqrt(C(1,:).^2+C(3,:).^2);

    % Linear terms of the equations of motion
    c = tf;                                                                         % Normalizing factor
    f = -[c.^2.*mu.*C(1,:)./r.^3; zeros(1,size(C,2)); c.^2.*mu.*C(3,:)./r.^3];      % Acceleration vector
    a = [C(7,:)-C(1,:).*C(5,:).^2; C(1,:).*C(8,:)+2*C(4,:).*C(5,:); C(9,:)];        % Inertial acceleration field
    dv = [C(4,:); C(1,:).*C(5,:); C(6,:)];                                          % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = a-f;
end
