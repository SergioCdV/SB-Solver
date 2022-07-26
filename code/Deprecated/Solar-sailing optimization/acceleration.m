%% Project: 
% Date: 30/01/22

%% Acceleration %%
% Function to compute the acceleration vector norm from cylindrical coordinates

% Inputs: - scalar mu, the gravitational parameter of the system 
%         - scalar r0, the characteristic or dimensionalising distance of
%         - scalar tf, the approximated time of flight of the mission
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - vector n, with the degrees of approximation of each position
%           coordinate

% Outputs: - scalar a, the acceleration vector norm

function [a] = acceleration(mu, r0, tf, P, B, n)
    % Extract the vector components of the state
    C = evaluate_state(P, B, n);
    [rho, v, gamma] = extract_coordinates(C, r0, tf);
    
    % Heliocentric position vector
    r = sqrt(rho(1,:).^2 + rho(3,:).^2);
    
    % Equations of motion
    arho = gamma(1,:) + mu.*rho(1,:)/r.^3 - rho(1,:).*v(2,:).^2;
    atheta = rho(1,:).*gamma(2,:) + 2.*v(1,:).*v(2,:);
    az = gamma(3,:) + mu.*rho(3,:)./r.^3;
    
    % Magnitude of the acceleration
    a = sqrt(arho.^2+ atheta.^2 + az.^2);
end