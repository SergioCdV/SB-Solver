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
    % State variables 
    a = C(1,:);             % Semimajor axis
    e = C(2,:);             % Eccentricity
    n = C(3,:);             % Mean motion
    v = C(4,:);             % True anomaly

    % Compute the atmospheric density
    rho = ones(size(a));

    % Compute the control vector (ballistic coefficient)
    vrel = n.^2..*a.^2./(1-e.^2).*(ones(size(e))+e.^2+2*e.*cos(v));
    E = sqrt(ones(size(e))+e.^2+2*e.*cos(v))./(n.*sqrt(ones(size(e))-e.^2));
    u = -C(5,:)./(rho.*vrel.^2.*E);
end