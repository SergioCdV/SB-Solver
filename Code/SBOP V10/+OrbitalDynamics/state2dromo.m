%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/04/22

%% Dromo elements to COE %%
% This file contains the function to change from dromo orbital elements to COE

% Inputs: - vector s, containing the inertial position and velocity vector. 

% Ouputs: - vector d, containing the state vector in the DROMO format.

function [d] = state2dromo(s)
    h = cross(s(1:3,1), s(4:6,1));                      % Angular momentum vector 
    e = -s(1:3,1)/norm(s(1:3,1))-cross(h, s(4:6,1));    % Eccentricity vector
    k = h/norm(h);                                      % Inertial k unit vector
    i = s(1:3,1)/norm(s(1:3,1));                        % Inertial i unit vector

    if (norm(e) >= 0 && norm(e) < 1)
        u_1 = e/norm(e);                                % Perifocal departure frame
        u_2 = cross(k, u_1);                            % Perifocal departure frame
        j = cross(k,i);                                 % Inertial j unit vector

        Q = [u_1, u_2, k].';                            % Initial rotation matrix of the ideal frame
        Omega = atan2(Q(3,1),-Q(3,2));                  % RAAN
        omega = atan2(Q(1,3),Q(2,3));                   % Argument of perigee
        inc = acos(Q(3,3));                             % Inclination of the orbit

        diff = (Omega-omega)/2;
        plus = (Omega+omega)/2;

        % Attitude quaternion
        q(1,1) = sin(inc/2) * cos(diff);
        q(2,1) = sin(inc/2) * sin(diff);
        q(3,1) = cos(inc/2) * sin(plus);
        q(4,1) = cos(inc/2) * cos(plus);
    
        sin_sigma = dot(-j, u_1);                       % Sine of the initial perifocal angle
        cos_sigma = dot(i, u_1);                        % Cosine of the initial perifocal angle
        sigma = atan2(sin_sigma, cos_sigma);            % Initial perifocal angle
        zeta = [norm(e) 0 1/norm(h)];                   % Dynamic variables
    
        % Final orbital elements
        d = [zeta.'; q; sigma];
    end
end