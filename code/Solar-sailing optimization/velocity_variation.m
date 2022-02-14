%% Project: 
% Date: 31/01/22

%% Extract coordinates %%
% Function to compute the total DV needed for a given trajectory, used as
% an optimisation objective function

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - scalar r0, the characteristic or dimensionalising distance of
%           the mission
%         - scalar T, the characteristic or dimensionalising time of
%           the mission
%         - vector x, the vector of optimization variables
%         - array B, the polynomial basis in use in the approximation
%         - vector n, the degree of approximation to each position
%           coordinate

% Outputs: - scalar DV, the total velocity change

function [DV] = velocity_variation(mu, r0, T, tau, x, B, n)
    % Extract the variables of the problem 
    P = reshape(x(1:end-1), [3, max(n)+1]);

    % Compute the magnitude of necessary propulsive acceleration in flight
    a = acceleration(mu, P, B, n, r0, T);

    % Trapezoidal approximation of the DV integral
    DV = trapz(tau, a)*(r0/T);
end
