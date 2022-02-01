%% Project: 
% Date: 31/01/22

%% Extract coordinates %%
% Function to compute the total DV needed for a given trajectory, used as
% an optimisation objective function

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - scalar r0, the characteristic or dimensionalising distance of
%         - scalar tf, the approximated time of flight of the mission
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation

% Outputs: - scalar DV, the total velocity change

function [DV] = velocity_variation(mu, r0, tau, tf, P, B)
    % Compute magnitude of necessary propulsive acceleration in flight
    a = acceleration(mu, r0, tf, P, B);

    % Trapezoidal approximation of the DV integral
    DV = trapz(tau, a);
end





