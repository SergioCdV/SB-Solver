%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the solution method
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
final_coe = [1.05*r0 5e-3 deg2rad(15) deg2rad(1) deg2rad(10)]; 
thetaf = deg2rad(270);
final_coe = [final_coe thetaf]; 

% Spacecraft propulsion parameters 
T = 0.05e-3;     % Maximum acceleration 

% Optimization 
method = 'Shape-based';

%% Results
% Simple solution    
tic
[t, S, dV] = MPC_wrapper(system, initial_coe, final_coe, T, method);
toc 
