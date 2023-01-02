%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Setup of the solution method
animations = 0;                         % Set to 1 to generate the gif
time_distribution = 'Chebyshev';        % Distribution of time intervals
basis = 'Chebyshev';                    % Polynomial basis to be use
n = 9;                         % Polynomial order in the state vector expansion
m = 100;                                % Number of sampling points
cost_function = 'Minimum fuel';         % Cost function to be minimized

% System data 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

system.mu = mu; 
system.distance = r0; 
system.time = t0; 

% Earth's orbital elements
initial_coe = [r0 1e-3 0 deg2rad(0) deg2rad(0)]; 
theta0 = deg2rad(0);
initial_coe = [initial_coe theta0]; 

% Mars' orbital elements 
final_coe = [1.05*r0 1e-3 deg2rad(5) deg2rad(50) deg2rad(0)]; 
thetaf = deg2rad(100);
final_coe = [final_coe thetaf]; 

% Spacecraft parameters 
T = 0.5e-3;              % Maximum acceleration 
TOF = 1*365*3600*24;     % Desired TOF for the time-fixed problem

% Setup 
setup.order = n; 
setup.basis = basis;
setup.grid = time_distribution; 
setup.nodes = m; 
setup.cost_function = cost_function;
setup.FreeTime = true;

setup.resultsFlag = true; 
setup.animations = false; 

%% Optimization
% Simple solution    
tic
[C, dV, u, tf, tfapp, tau, thetaf, exitflag, output] = sb_solver(system, initial_coe, final_coe, TOF, T, setup);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, tf, tfapp, tau, thetaf, exitflag, output] = sb_solver(system, initial_coe, final_coe, TOF, T, setup);
    time(i) = toc;
end

time = mean(time);
