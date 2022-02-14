%% Project: 
% Date: 31/01/22

%% Initial data %%
% Function to generate the boundary conditions of interest 

% Inputs: - scalar r0, the characteristic distance of the mission
%         - scalar T, the characteristic time
%         - scalar cases, the case study to be studied

% Outputs: - vector initial, the initial state vector 
%          - vector final, the final state vector

function  [initial, final] = initial_data(r0, T, cases)
    switch (cases)
        case 1
            [initial, final] = case_1(r0,T);
        case 2 
            [initial, final] = case_2(r0,T);
        otherwise 
            error('No valid mission was selected');
    end
end

%% Auxiliary functions
% Test case 1 data. Source: https://ssd.jpl.nasa.gov/horizons.cgi
function [initial, final] = case_1(r0,T)
    % Preallocation for correct shaping 
    initial = zeros(11,1);                          % Initial boundary conditions
    final = zeros(11,1);                            % Final boundary conditions

    % Spacecraft initial cylindrical position vector
    initial(1) = 0.99168*r0;                        % In-plane radial coordinate [m]
    initial(2) = 2.8169;                            % Longitude [rad]
    initial(3) = 1.085074459426786E-04*r0;          % Altitude [m]
    
    % Spacecraft initial cylindrical velocity vector
    initial(4) = 0;                                 % In-plane radial coordinate velocity [m/s]
    initial(5) = 29.87e3/initial(1);                % Longitude velocity [rad/s]
    initial(6) = 3e3/initial(1);                    % Altitude velocity [m/s]
    
    % Spacecraft initial orbital elements (semimajor axis and inclination)
    initial(10) = r0;                               % Semimajor axis [m]
    initial(11) = 1e-2;                             % Inclination [rad]
    
    % Final spacecraft in-plane Cartesian coordinates
    xf = 0.47022191181092;
    yf = 1.44869360058;
    
    % Final spacecraft cylindrical position vector
    final(1) = 1.523213981040*r0;                   % In-plane radial coordinate [m]
    final(2) = 4*pi + atan2(yf,xf);                 % Longitude [rad] 
    final(3) = 0.0188*r0;                           % Altitude [m]
    
    % Final spacecraft cylindrical velocity vector
    final(4) = 0;                                   % In-plane radial coordinate velocity [m/s]
    final(5) = 24.07e3/final(1);                    % Longitude velocity [rad/s]
    final(6) = 1e-5;                                % Altitude velocity [m/s]
    
    % Spacecraft final orbital elements
    final(10) = final(1);                           % Semimajor axis [m]
    final(11) = deg2rad(1.85) - initial(11);        % Inclination [rad]
end

% Test case 2 data. Source: https://ssd.jpl.nasa.gov/horizons.cgi
function [initial, final] = case_2(r0,T)
    % Preallocation for correct shaping 
    initial = zeros(11,1);                          % Initial boundary conditions
    final = zeros(11,1);                            % Final boundary conditions

    % Spacecraft initial cylindrical position vector
    initial(1) = 0.99168*r0;                        % In-plane radial coordinate [m]
    initial(2) = 2.8169;                            % Longitude [rad]
    initial(3) = 1.085074459426786E-04*r0;          % Altitude [m]
    
    % Spacecraft initial cylindrical velocity vector
    initial(4) = 0;                                 % In-plane radial coordinate velocity [m/s]
    initial(5) = 29.87e3/initial(1);                % Longitude velocity [rad/s]
    initial(6) = 3e3/initial(1);                    % Altitude velocity [m/s]
    
    % Spacecraft initial orbital elements (semimajor axis and inclination)
    initial(10) = r0;                               % Semimajor axis [m]
    initial(11) = 1e-2;                             % Inclination [rad]
    
    % Final spacecraft in-plane Cartesian coordinates
    xf = 0.47022191181092;
    yf = 1.44869360058;
    
    % Final spacecraft cylindrical position vector
    final(1) = 1.523213981040*r0;                   % In-plane radial coordinate [m]
    final(2) = 4*pi + atan2(yf,xf);                 % Longitude [rad] 
    final(3) = 0.0188*r0;                           % Altitude [m]
    
    % Final spacecraft cylindrical velocity vector
    final(4) = 0;                                   % In-plane radial coordinate velocity [m/s]
    final(5) = 24.07e3/final.pos(1);                % Longitude velocity [rad/s]
    final(6) = 1e-5;                                % Altitude velocity [m/s]
    
    % Spacecraft final orbital elements
    final(10) = final(1);                           % Semimajor axis [m]
    final(11) = deg2rad(1.85) - initial(11);        % Inclination [rad]
end

