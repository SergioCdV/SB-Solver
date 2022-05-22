%% Project: 
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the collocation method
time_distribution = 'Chebyshev';         % Distribution of time intervals
basis = 'Chebyshev';                     % Polynomial basis to be use
n = 12;                                  % Order of Bezier curve functions for each coordinate
m = 100;                                  % Number of sampling points

% System data 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

system.mu = mu; 
system.distance = r0; 
system.time = t0; 

% Spacecraft propulsion parameters 
T = 0.5e-3;     % Maximum acceleration 

% Initial input revolutions 
K = 0;

% Setup 
setup.resultsFlag = false; 
setup.animations = false; 

% Earth's orbital elements
initial_coe = [r0 1e-4 0 deg2rad(0) 0]; 
theta0 = deg2rad(0);
initial_coe = [initial_coe theta0]; 
Q0 = [cos(theta0) sin(theta0) 0; -sin(theta0) cos(theta0) 0; 0 0 1];

%% Monte-carlo analysis 
cases = 100;                       % Number of Monte-carlo cases to consider

% Preallocation for speed
Def_factor = zeros(2, cases);
V = zeros(1, cases);
Tf = zeros(1, cases);
S = zeros(1, cases);

for i = 1:cases
    % Orbital elements' statistical distribution 
    a = 0.9+0.2*rand;             % Semimajor axis distribution 
    e = 1e-4+1e-2*rand;           % Eccentricity distribution 
    RAAN = deg2rad(20*rand);      % RAAN distribution
    I = deg2rad(5*rand);          % Inclination distribution
    omega = deg2rad(360*rand);    % AoP distribution
    thetaf = deg2rad(360*rand);   % True anomaly distribution 

    % Final orbital elements
    final_coe = [a*r0 e RAAN I omega thetaf];
    Qf = euler_matrix(final_coe); % Final Euler matrix

    % Optimisation
    [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, time_distribution, basis, n, setup);
    %[sol, C, dV, u, tf, tfapp, tau, exitflag, output] = ga_wrapper(system, initial_coe, final_coe, K, T, setup);

    % Results
    da = initial_coe(1)-final_coe(1);
    db = initial_coe(1)*sqrt(1-initial_coe(2)^2)-final_coe(1)*sqrt(1-final_coe(2)^2);
    Def_factor(1,i) = norm([da db]) / initial_coe(1);
    Def_factor(2,i) = acos((1/2)*(trace(Q0*Qf.')-1));

    V(i) = dV;                    % Final dV cost
    Tf(i) = tf;                   % Final time of flight

    % Optimization success
    S(i) = (exitflag == 1) || (exitflag == 2) ||(exitflag == 3);       
end

%% Results
% 3D plot 
figure 
scatter3(Def_factor(1,:), Def_factor(2,:), V, 'filled'); 
grid on; 
xlabel('Orbit size difference $\Delta a$')
ylabel('Orbit orientation difference $\delta$')
zlabel('Transfer cost $\Delta V$ [m/s]')

figure 
scatter3(Def_factor(1,:), Def_factor(2,:), Tf, 'filled'); 
grid on; 
xlabel('Orbit size difference $\Delta a$')
ylabel('Orbit orientation difference $\delta$')
zlabel('Transfer time $t_f$ [s]');

% Porckchop plot 
figure 
% Sort results 
[Def_factor(1,:), index] = sort(Def_factor(1,:));
V = V(index); 
Tf = Tf(index);
subplot(1,2,1)
plot(Def_factor(1,:), V);
grid on; 
xlabel('Orbit size difference $\Delta a$')
ylabel('Transfer cost $\Delta V$ [m/s]')

[Def_factor(2,:), index] = sort(Def_factor(2,:));
V = V(index); 
Tf = Tf(index);
subplot(1,2,2)
plot(Def_factor(2,:), V);
grid on;
xlabel('Orbit orientation difference $\delta$')
ylabel('Transfer cost $\Delta V$ [m/s]')
sgtitle('Transfer porckchop plot')

figure
% Sort results 
[Def_factor(1,:), index] = sort(Def_factor(1,:));
V = V(index); 
Tf = Tf(index);
subplot(1,2,1)
plot(Def_factor(1,:), Tf);
grid on; 
xlabel('Orbit size difference $\Delta a$')
ylabel('Transfer time $t_f$ [s]')

[Def_factor(2,:), index] = sort(Def_factor(2,:));
V = V(index); 
Tf = Tf(index);
subplot(1,2,2)
plot(Def_factor(2,:), Tf);
grid on; 
xlabel('Orbit orientation difference $\delta$')
ylabel('Transfer time $t_f$ [s]')
sgtitle('Transfer time plot')