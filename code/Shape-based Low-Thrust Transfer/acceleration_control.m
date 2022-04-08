%% Project: 
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the 9xm state vector 
%         - scalar tf, the final time of flight
% Outputs: - vector u, the nondimensional 3xm control vector

function [u] = acceleration_control(mu,C,tf)
    % Compute the control vector
    r = sqrt(C(1,:).^2+C(3,:).^2);
    u = [C(7,:)+(-C(1,:).*C(5,:).^2+tf^2*mu.*C(1,:)/r.^3); ...
         C(1,:).*C(8,:)+(2*C(4,:).*C(5,:)); ... 
         C(9,:)+tf^2*(mu.*C(3,:)/r.^3)];
end