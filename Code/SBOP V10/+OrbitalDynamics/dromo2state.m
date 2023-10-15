%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/04/22

%% Dromo elements to Cartesian state vector %%
% This file contains the function to change from dromo orbital elements to Cartesian state vector

% Inputs: - vector elements, containing the DROMO elements. 

% Ouputs: - vector s, containing the state vector in COE format.

function [s] = dromo2state(elements)
    S = 1+elements(1).*cos(elements(8))+elements(2).*sin(elements(8));

    % Osculating motion plane position and velocity
    r = [cos(elements(8)); sin(elements(8)); 0] / (elements(3)^2*S);
    v = [-sin(elements(8)) - elements(2); cos(elements(8)) + elements(1); 0] * elements(3);

    % Inertial position and velocity 
    r = QuaternionAlgebra.RotateVector(QuaternionAlgebra.quaternion_inverse(elements(4:7)), r);
    v = QuaternionAlgebra.RotateVector(QuaternionAlgebra.quaternion_inverse(elements(4:7)), v);

    % Final state vector
    s = [r; v];
end