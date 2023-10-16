%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/04/22

%% COE to DROMO elements %%
% This file contains the function to change from COE to dromo orbital elements

% Inputs: - vector elements, containing the COE elements

% Ouputs: - vector s, containing the state vector in DROMO format

function [s] = coe2dromo(mu, elements)
    % Dynamic variables 
    zeta = [elements(2); 0; 1/sqrt(mu * elements(1) * (1-elements(2)^2))];
    
    % Attitude quaternion
    Omega = elements(3);
    inc = elements(4);
    omega = elements(5);

    diff = (Omega-omega)/2;
    plus = (Omega+omega)/2;
    q(1,1) = sin(inc/2) * cos(diff);
    q(2,1) = sin(inc/2) * sin(diff);
    q(3,1) = cos(inc/2) * sin(plus);
    q(4,1) = cos(inc/2) * cos(plus);

    beta = atan2(q(2,1), q(1,1));
    
    % Final DROMO elements
    s = [zeta; q; elements(6)+beta];
end