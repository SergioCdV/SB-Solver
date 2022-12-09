%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Classical orbital elements to equinoctial %% 
% Function to transform classical orbital elements into Cartesian or
% viceversa

% Inputs: - vector s, the orbital state vector to be transformed 
%         - boolean direction, 0 for the COE to equinoctional
%           transformation and 1 for the viceversa case

% Outputs: - vector S, the transformed orbital state vector

function [S] = coe2equinoctial(s, direction)
    % Sanity check on the s dimensions 
    if (size(s,1) ~= 6)
        lastwarn('Input orbital element set is not valid')
    end

    % Switch directions 
    if (direction)
        S(1) = s(1)*(1-s(2)^2);             % Semiparameter
        S(2) = s(2)*cos(s(5)+s(3));
        S(3) = s(2)*sin(s(5)+s(3));
        S(4) = tan(s(4)/2)*cos(s(3));
        S(5) = tan(s(4)/2)*sin(s(3));
        S(6) = s(6)+s(5)+s(3);
    
    else
        S(1) = s(1)/(1-s(2)^2-s(3)^2);                          % Semimajor axis
        S(2) = sqrt(s(2)^2+s(3)^2);                             % Eccentricity
        S(3) = atan2(s(5),s(4));                                % RAAN
        S(4) = atan2(2*sqrt(s(4)^2+s(5)^2), 1-s(4)^2-s(5)^2);   % Inclination
        S(5) = atan2(s(3)*s(4)-s(2)*s(5),s(2)*s(4)+s(3)*s(5));  % Argument of periapsis
        S(6) = s(6)-(S(5)+S(3));                                % True anomaly
    end  
end