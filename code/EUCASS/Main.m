%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the solution method
time_distribution = 'Chebyshev';      % Distribution of time intervals
basis = 'Chebyshev';                    % Polynomial basis to be use
n = [20 20 20];                          % Polynomial order in the state vector expansion
m = 600;                                 % Number of sampling points

% System data 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

system.mu = mu; 
system.distance = r0; 
system.time = t0; 

% Earth's orbital elements
initial_coe = [r0 1e-3 0 deg2rad(0) deg2rad(0)]; 
theta0 = deg2rad(95);
initial_coe = [initial_coe theta0]; 

% Mars' orbital elements 
final_coe = [1.05*r0 1e-4 deg2rad(5) deg2rad(1) deg2rad(30)]; 
thetaf = deg2rad(270);
final_coe = [final_coe thetaf]; 

% Spacecraft propulsion parameters 
T = 0.05e-3;     % Maximum acceleration 

% Initial input revolutions 
K = 1;

% Setup 
setup.resultsFlag = true; 
setup.animations = false; 

%% Results
% Simple solution    
tic
[C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, time_distribution, basis, n, setup);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, time_distribution, basis, n, setup);
    time(i) = toc;
end

time = mean(time);
