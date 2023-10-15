%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/04/22

%% Dromo elements to COE %%
% This file contains the function to change from dromo orbital elements to COE

% Inputs: - vector elements, containing the DROMO elements. 

% Ouputs: - vector s, containing the state vector in COE format.

function [s] = dromo2coe(elements)
    e = sqrt(elements(1)^2 + elements(2)^2);    % Eccentricity
    a = -1/((e^2-1) * elements(3)^2);           % Semimajor axis 

    beta = atan2(elements(2), elements(1));
    i = 2 * acos(sqrt(elements(6)^2 + elements(7)^2));   % Orbital inclination
    p = 2 * atan2(elements(6), elements(7));             % RAAN + AoP
    r = 2 * atan2(elements(5), elements(4));             % RAAN - AoP
    
    Omega = p + r;                                       % RAAN
    omega_tilde = Omega - r;
    omega = beta + omega_tilde;                          % AoP
    
    % Final orbital elements
    s = [a e i Omega omega elements(8)-beta];
end