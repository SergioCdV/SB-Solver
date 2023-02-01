%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Setup of the solution method
time_distribution = 'Legendre';        % Distribution of time intervals
basis = 'Legendre';                    % Polynomial basis to be use
n = [10 12 12];                        % Polynomial order in the state vector expansion
m = 60;                                % Number of sampling points

% Earth's orbital elements
initial_coe = [r0 1e-3 0 deg2rad(0) deg2rad(0)]; 
theta0 = deg2rad(0);
initial_coe = [initial_coe theta0]; 

% Mars' orbital elements 
final_coe = [7.10*r0 1e-2 deg2rad(0) deg2rad(20) deg2rad(5)]; 
thetaf = deg2rad(20);
final_coe = [final_coe thetaf]; 

% Spacecraft parameters 
T = 0.5e-4;              % Maximum acceleration 
TOF = 1*365*3600*24;     % Desired TOF for the time-fixed problem

% Setup 
setup.order = n; 
setup.basis = basis;
setup.grid = time_distribution; 
setup.nodes = m; 
setup.cost_function = cost_function;
setup.FreeTime = true;

setup.resultsFlag = true; 

%% Optimization
% Simple solution    
tic
[C, cost, u, tf, tfapp, tau, exitflag, output] = sb_solver(initial_state, final_state, TOF, T, setup);
toc 
