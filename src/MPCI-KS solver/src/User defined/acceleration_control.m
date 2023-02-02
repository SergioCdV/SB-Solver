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
    % Linear terms of the equations of motion
    a = C(11:14,:);                                                                 % Inertial acceleration field
    f = -tf^2*(mu.*C(5,:)/4).*C(1:4,:);                                             % Acceleration vector
    dv = C(6:9,:);                                                                  % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = a-f;
    for i = 1:size(C,2)
        L = KS_matrix(C(1:4,i));
        u(:,i) = L*u(:,i);
    end

    u = 2*u./dot(C(1:4,:),C(1:4,:),1).^2;
end
