%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Setup of the solution method
animations = 0;                         % Set to 1 to generate the gif
time_distribution = 'Bernstein';        % Distribution of time intervals
basis = 'Bernstein';                    % Polynomial basis to be use
n = [13 13 13];                         % Polynomial order in the state vector expansion
m = 500;                                % Number of sampling points
cost_function = 'Minimum fuel';         % Cost function to be minimized

% System data 
M1 = 1.9891e+30;                        % Mass of the Sun [kg] 
M2 = 5.97219e+24;                       % Mass of the Earth [kg] 
r0 = 149597870700;                      % 1 AU [m]
t0 = 86400*365;                         % Fundamental time unit

M1 = 5.97219e+24;                       % Mass of the Earth [kg] 
M2 = 7.34767309e+22;                    % Mass of the Moon [kg] 
r0 = 3.850e8;                           % Earth-Moon distance [m]
t0 = 86400*28;                          % Fundamental time unit

system.M1 = M1; 
system.M2 = M2; 
system.distance = r0; 
system.time = t0; 

% Initial ECI position and velocity
initial = [37e6 0 0 0 sqrt(3.986e14/(37e6)) 0]; 

% Final ECI position and velocity
final = [0.98*r0 0 0 -sqrt(3.986e14/(0.98*r0)) 0 0]; 

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
[C, dV, u, tf, tfapp, tau, exitflag, output] = sb_solver(system, initial, final, TOF, T, setup);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, tf, tfapp, tau, exitflag, output] = sb_solver(system, initial, final, TOF, T, setup);
    time(i) = toc;
end

time = mean(time);
