%% Project: 
% Date: 31/01/22

%% Extract coordinates %%
% Function to compute the total DV needed for a given trajectory, used as
% an optimisation objective function

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - scalar r0, the characteristic or dimensionalising distance of
%         - vector x, the vector of optimization variables
%         - array B, the polynomial basis in use in the approximation
%         - vector n, the degree of approximation to each position
%           coordinate

% Outputs: - scalar DV, the total velocity change

function [DV] = velocity_variation(mu, r0, tau, x, B, n)
    % Extract the variables of the problem 
    P = reshape(x(1:end-1), [3, max(n)+1]);
    tf = x(end);

    % Compute the magnitude of necessary propulsive acceleration in flight
    a = acceleration(mu, r0, tf, P, B, n);

    % Trapezoidal approximation of the DV integral
    DV = trapz(tau, a)*tf;
end





