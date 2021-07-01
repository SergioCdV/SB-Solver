function [initial, final] = initial_coordinates(r0)

% INITIAL_COORDINATES function defines the initial and final position and
% velocity of the solar sail vehicle, with trial case defined.

% source: https://ssd.jpl.nasa.gov/horizons.cgi


% initial coordinates of spacecraft
initial.pos(1) = 0.99168*r0;         % [m]
initial.pos(2) = 2.8169;          % [rad]
initial.pos(3) = 0; %1.085074459426786E-04*r0;          % [m]

% initial velocity
initial.vel(1) = 0;                        % [m/s]
initial.vel(2) = 29.87e3/initial.pos(1);      % [rad/s]
initial.vel(3) = 0;                        % [m/s]

% orbital elements
initial.a = r0;
initial.i = 0; %1e-2;

% final coordinates of spacecraft
y0 = 0.31635;
x0 = 0.9398690;

yf = 1.44869;
xf = 0.470221;

final.pos(1) = 1.5232*r0;         % [m]
final.pos(2) = 3*pi + atan(y0/x0) + atan(yf/xf)  + 2.8169;       % one and a half full orbits, plus the two extra sections [rad] 
final.pos(3) = 0.0188*r0;        % [m]


% final velocity
final.vel(1) = 0;                             % [m/s]
final.vel(2) = 24.07e3/final.pos(1);      % [rad/s]
final.vel(3) = 1e-5;                             % [m/s]

% orbital elements
final.a = final.pos(1);
final.i = deg2rad(1.85) - initial.i;


end