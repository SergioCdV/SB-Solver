function [initial, final] = initial_coordinates2(r0)

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

% jul-23
xf = 5.586904658990008E-01;
yf = 1.407225129071290E+00;
zf = 1.591100056503527E-02; 

% jul-23
xf = 0.47022191181092;
yf = 1.44869360058;
zf = 0.0188;

% jul-31
xf = 3.667890234308765E-01;
yf = 1.488822261886027E+00;
zf = 2.232422514156584E-02;

% % aug-07
% xf = 2.746948002856574E-01;
% yf = 1.517538314367233E+00;
% zf = 2.518309270020771E-02;



final.pos(1) = sqrt(yf^2+xf^2)*r0;         % [m]
final.pos(2) = 4*pi + atan2(yf,xf);       % one and a half full orbits, plus the two extra sections [rad] 
final.pos(3) = zf*r0;        % [m]


% final velocity
final.vel(1) = 0;                             % [m/s]
final.vel(2) = 24.07e3/final.pos(1);      % [rad/s]
final.vel(3) = 0;                             % [m/s]

% orbital elements
final.a = final.pos(1);
final.i = deg2rad(1.85) - initial.i;


end