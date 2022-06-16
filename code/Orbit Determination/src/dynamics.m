%% Project: 
% Date: 16/06/22

%% Dynamics %%
% Function to compute the dynamics vector field of the spacecraft

% Inputs: - scalar mu, the gravitational parameter of the orbited central
%           body
%         - scalar J2, the second zonal harmonic of the central body 
%         - scalar Re, the central body reference radius
%         - vector s, the Cartesian state vector of the spacecraft 
%         - vector u, the acceleration control law of the spacecraft

% Outputs: - vector field ds, describing the dynamics under consideration

function [ds] = dynamics(mu, J2, Re, t, s, u)
    % Disassemble the state vector
    r = s(1:3);         % ECI position vector
    v = s(4:6);         % ECI velocity vector

    % Compute the Newtonian term 
    an = -(mu/norm(r)^3)*r;

    % Compute the J2 perturbation term 
    ap = mu*J2*(Re^2/norm(r)^7)*[r(1)*(6*r(3)^2-(3/2)*(r(1)^2+r(2)^2)); ...
                                 r(2)*(6*r(3)^2-(3/2)*(r(1)^2+r(2)^2)); ...
                                 r(3)*(3*r(3)^2-(9/2)*(r(1)^2+r(2)^2))];
    ap = zeros(3,1);

    % Total acceleration applied to the body
    a = an+ap+u; 

    % Assemble the state vector derivative 
    ds = [v; a];
end