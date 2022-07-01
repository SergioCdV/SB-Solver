%% Project: Shape-based attitude planning %%
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - matrix I, the inertia dyadic of the system
%         - array C, the 9xm state vector 
%         - scalar tf, the final time of flight

% Outputs: - vector u, the nondimensional 3xm control vector

function [u] = acceleration_control(I, C, tf)
    % Compute the control vector as a residual of the dynamics
    c = tf;             % Normalizing factor

    % Pre-allocation 
    u = zeros(3,size(C,2)); 

    % Compute the control vector as a dynamics residual
    for i = 1:size(C,2)
        % Dynamics
        omega = 2*quaternion_product(C(5:8,i), quaternion_inverse(C(1:4,i)));
        Domega = 2*quaternion_product(C(9:12,i),quaternion_inverse(C(1:4,i)))+[0; cross(omega(2:4),I*omega(2:4))];
        u(:,i) = Domega(2:4);
    end
end
