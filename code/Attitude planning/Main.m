%% Project: 
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the solution method
time_distribution = 'Linear';      % Distribution of time intervals
basis = 'Bernstein';               % Polynomial basis to be use
n = 9;                             % Order of Bezier curve functions for each coordinate
m = 60;                            % Number of sampling points

% System data 
I = diag([1 1 1]);     % Inertia dyadic of the system  

system.Inertia = I; 
system.V = [0;1;0];
system.Prohibited = [sqrt(2)/2;sqrt(2)/2;0];
system.Tol = 1e-1; 

% Initial conditions (Euler angles + angular velocity)
initial_bc = [deg2rad(90) 0 0 zeros(1,3)];  

% Final conditions (Euler angles + angular velocity)
final_bc = [deg2rad(45) deg2rad(90) deg2rad(300) zeros(1,3)]; 

% Maximum applicable torque
T = 0.05e-3;    

% Setup 
setup.resultsFlag = true; 
setup.animations = false; 

%% Results
% Simple solution    
tic
[C, domega, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_bc, final_bc, T, m, time_distribution, basis, n, setup);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_bc, final_bc, K, T, m, time_distribution, basis, n, setup);
    time(i) = toc;
end

time = mean(time);