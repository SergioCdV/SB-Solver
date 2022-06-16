%% Project: 
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
initial_coe = [r0 1e-3 0 deg2rad(0) 0]; 
theta0 = deg2rad(95);
initial_coe = [initial_coe theta0]; 

% Mars' orbital elements 
final_coe = [1.05*r0 5e-3 deg2rad(15) deg2rad(1) deg2rad(10)]; 
thetaf = deg2rad(270);
final_coe = [final_coe thetaf]; 

% Spacecraft propulsion parameters 
T = 0.05e-3;     % Maximum acceleration 

% Initial input revolutions 
K = 1;

% Setup 
setup.resultsFlag = false; 
setup.animations = false; 

%% Optimization 
[sol, fval] = ga_wrapper(system, initial_coe, final_coe, K, T, setup);

%% Results
dV = zeros(1,size(sol.NumPoints,1));
tf  = dV;
time_av = tf;

for i = 1:size(sol.NumPoints,2)
    m = sol.NumPoints(i); 
    n = sol.Order(i); 
    basis = sol.Basis{i}; 
    time_distribution = sol.Points{i}; 

    setup.resultsFlag = true; 
    tic
    [C, dV(i), u, tf(i), tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, time_distribution, basis, n, setup);
    toc 

    % Average results 
    iter = 0; 
    time = zeros(1,iter);
    setup.resultsFlag = false; 
    for j = 1:iter
        tic 
        [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, time_distribution, basis, n, setup);
        time(j) = toc;
    end
    
    time_av(i) = mean(time);
end

%% Plots 
figure 
plot(fval(:,1), fval(:,2), 'o')
xlabel('Normalized transfer cost $\Delta V$')
ylabel('Normalized time of flight $t_f$')
grid on;

figure 
plot(sol.NumPoints, sol.Order, 'o')
xlabel('Number of sampling nodes $m$')
ylabel('Degree of expansion $n$')
grid on;