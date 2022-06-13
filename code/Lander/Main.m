%% Project: 
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the collocation method
time_distribution = 'Linear';            % Distribution of time intervals
basis = 'Chebyshev';                     % Polynomial basis to be use
n = 9;                                   % Order of Bezier curve functions for each coordinate
m = 20;                                 % Number of sampling points

% System data
a = [490 490 490];                      % Semimajor axes of the ellipsoid
r0 = max(a);                            % 1 AU [m]
mu = 2*pi*6.67e-11*7.33e10;             % Gravitational parameter of Bennu [m^3 s^−2] 
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

system.mu = mu; 
system.distance = r0; 
system.time = t0; 
system.ellipsoid  = a; 

% Initial landing conditions
initial = [0 0 2*a(3) 0 0 -10];  

% Final landing conditions
final = [0 0 a(3) zeros(1,3)]; 

% Spacecraft propulsion parameters 
T = 1e3;     % Maximum acceleration 

% Setup 
setup.resultsFlag = true; 
setup.animations = false; 

%% Optimization 
[C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial, final, T, m, time_distribution, basis, n, setup);
