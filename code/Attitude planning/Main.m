%% Project: Shape-based attitude planning %%
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                        % Set to 1 to generate the gif

%% Setup of the solution method
time_distribution = 'Chebyshev';      % Distribution of time intervals
basis = 'Chebyshev';                  % Polynomial basis to be use
n = 5;                               % Order of approximation functions for each coordinate
m = 100;                              % Number of sampling points
cost = 'Minimum time';              % Cost function to be minimized 
dynamics = 'Euler';                   % Formulation of the dynamics vectorfield to be used

% System data 
I = diag([1 2 3]);                    % Inertia dyadic of the system  

system.Inertia = I; 
system.V = [0;1;0];
system.Prohibited = [sqrt(2)/2; 0; sqrt(2)/2];
system.Tol = 1e-1; 

% Initial conditions (Euler angles + angular velocity)
initial_bc = [deg2rad(0) 0 0 zeros(1,3)];  

% Final conditions (Euler angles + angular velocity)
final_bc = [deg2rad(90) deg2rad(0) deg2rad(270) ones(1,3)]; 

% Maximum applicable torque
T = 1e-6;    

% Maneuver time 
tf = 600; 

% Setup 
setup.order = n;                        % Order in the approximation of the state vector
setup.formulation = dynamics;           % Formulation of the dynamics
setup.basis = basis;                    % Polynomial basis to be used 
setup.grid = time_distribution;         % Sampling grid to be used
setup.nodes = m;                        % Number of nodes in the grid
setup.cost_function = cost;             % Cost function to be minimized
setup.resultsFlag = true; 
setup.animations = false; 

%% Results
% Simple solution    
tic
[C, domega, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_bc, final_bc, tf, T, setup);
toc 

%%
% MPC solution
method = 'Shape-based';
setup.resultsFlag = false;
tic 
[t, S, U, dV] = MPC_wrapper(system, initial_bc, final_bc, T, method, setup);
toc

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_bc, final_bc, K, T, setup);
    time(i) = toc;
end

time = mean(time);
