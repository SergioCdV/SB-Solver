%% Project: 
% Date: 01/02/22

%% Orbital elements %%
% Function to compute the classical orbital elements (semimajor axis and
% inclination) of a given orbit

% Inputs: - scalar mu, the orbited body gravitational parameter
%         - vector s, the cylindrical state vector


% Outputs: - mu, the gravitational parameter of the central body 
%          - s, the state vector defining the orbit

function [a, i] = orbital_elements(mu, s)
    % Extract data and convert it to Cartesian coordinates
    S = cylindrical2cartesian(s, true);
    r = S(1:3);                                 % Spacecraft position vector
    v = S(4:6);                                 % Spacecraft velocity vector
    
    % Calculate the angular momentum vector
    h = cross(r,v);
    
    % Calculate the eccentricity vector
    e = cross(v,h)./mu - r./norm(r);
    
    % Calculate the inclination
    i = acos(h(3)/norm(h));
    
    % Calculate the semi-major axis
    a = norm(h)^2/mu/(1-norm(e)^2);

    % Compute the Euler matrix 
    Q = [e/norm(e); cross(h/norm(h), e/norm(e)); h/norm(h)];
end