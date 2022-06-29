%% Project: 
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the solution method
time_distribution = 'Regularized';      % Distribution of time intervals
basis = 'Bernstein';                    % Polynomial basis to be use
n = 10;                                 % Order of Bezier curve functions for each coordinate
m = 60;                                 % Number of sampling points

% Time distributions plot 
y = [zeros(1,m); 0.25*ones(1,m); 0.5*ones(1,m); ones(1,m)];
tau = [collocation_grid(m, 'Linear', ''); collocation_grid(m, 'Normal', ''); collocation_grid(m, 'Random', ''); collocation_grid(m, 'Chebyshev', '')];
figure 
ax = gca;
ax.YAxis.Visible = 'off';
hold on 
plot(tau(1,:), y(1,:), '-*');
plot(tau(2,:), y(2,:), '-*');
plot(tau(3,:), y(3,:), '-*');
plot(tau(4,:), y(4,:), '-*');
xlabel('$\tau_j$');
legend('Linear', 'Normal', 'Random', 'Chebyshev')

% a = [1011.44 935.65 1010.64 1009.13 1010.96 1010.22 999.70 1010.26 1136.49 1065.70]; 
% t = [587.90 594.57 588.78 613.45 861.06 587.88 599.72 587.86 589.27 588.96];
% 
% dVm = mean(a)+3*std(a);
% T = mean(t)+3+std(t); 

% System data 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

a = [1011.44 935.65 1010.64 1009.13 1010.96 1010.22 999.70 1010.26 1136.49 1065.70]/(r0/t0); 
t = 86400*[587.90 594.57 588.78 613.45 861.06 587.88 599.72 587.86 589.27 588.96]/t0;

system.mu = mu; 
system.distance = r0; 
system.time = t0; 

% Earth's orbital elements
initial_coe = [r0 1e-3 0 deg2rad(0) 0]; 
theta0 = deg2rad(100);
initial_coe = [initial_coe theta0]; 

% Mars' orbital elements 
final_coe = [1.05*r0 5e-3 deg2rad(5) deg2rad(1) deg2rad(5)]; 
thetaf = deg2rad(355);
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
