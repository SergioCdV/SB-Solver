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
    % Compute the auxiliary variables
    w = 1+C(2,:).*cos(C(6,:))+C(3,:).*cos(C(6,:));

    % Linear terms of the equations of motion
    f = [zeros(5,size(C,2)); sqrt(mu.*C(1,:)).*(w./C(1,:)).^2];                     % Acceleration vector
    a = C(7:12,:);                                                                  % Inertial acceleration field
    dv = C(7:12,:);                                                                 % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = zeros(3,size(C,2));

    for i = 1:size(C,2)
        B = control_input(mu, C(:,i));
        b = B([2 3 4],:);
        u(:,i) = b\(a([2 3 4],i)/tf-f([2 3 4],i));
    end
end


