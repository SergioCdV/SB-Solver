%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/04/22

%% Dromo elements to Cartesian state vector %%
% This file contains the function to change from dromo orbital elements to Cartesian state vector

% Inputs: - vector elements, containing the DROMO elements. 

% Ouputs: - vector s, containing the state vector in COE format.

function [s] = dromo32state(elements)
    l = 2*elements(1)-elements(1)^2-elements(2)^2;

    % Osculating motion plane position and velocity
    r = [elements(3)*l/elements(1); 0; 0];
    v = [elements(2); elements(1); 0] / (l*elements(3));

    % Inertial position and velocity 
    r = QuaternionAlgebra.RotateVector(QuaternionAlgebra.quaternion_inverse(elements(4:7)), r);
    v = QuaternionAlgebra.RotateVector(QuaternionAlgebra.quaternion_inverse(elements(4:7)), v);

    % Final state vector
    s = [r; v];
end