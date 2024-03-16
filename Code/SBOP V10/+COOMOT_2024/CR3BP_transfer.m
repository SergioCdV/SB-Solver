%% Project: SBOPT %%
% Date: 09/03/24

%% 3D low-thrust transfer in the CR3BP problem %% 
% This script provides a main interface to solve 3D low-thrust transfers in the CR3BP %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 15;                                % Polynomial order in the state vector expansion
m = 200;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

% System data 
mu = 1.2150581623433623E-2;     % Gravitational parameter of the Earth
Lem = 384400e3;                 % Mean distance from the Earth to the Moon
T0 = 28 * 86400;                % Characteristic time of the Earth-Moon system
n = 2*pi / T0;                  % Characteristic frequency of the system
Vc = Lem * n;                   % Characteristic velocity of the Earth-Moon system
gamma = Lem * n^2;              % Characteristic acceleration of the Earth-Moon system

% Initial and final clocks 
t0 = 0;                         % Initial clock 
tf = 2 * pi;                    % Final clock
    
% Initial and final conditions 
S0 = [0.9261 0 0.3616 0 -0.0544  0].';    % State vector of a vertical orbit                 
SF = [1.0406 0 0.1735 0 -0.0770 0].';     % State vector of a butterfly orbit; 

% Mission parameters 
T = 0.5e-3;                     % Maximum acceleration 
T = T/gamma;                    % Normalized acceleration

problem_params = [mu; t0; tf; T];

% Create the problem
OptProblem = Problems.CR3BP_transfer(S0, SF, L, StateDimension, ControlDimension, problem_params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

time = mean(time);

dV = dV * Vc;

%% Results 
% Configuration space coordinates
x = C(1,:);         % Synodic x coordinate
y = C(2,:);         % Synodic y coordinate
z = C(3,:);         % Synodic z coordinate
    
%% Plots
% Orbit representation
figure_orbits = figure;
view(3)
hold on
xlabel('$X$')
ylabel('$Y$')
zlabel('$Z$')
plot3(0,0,0, '*k');
plot3(x(1), y(1), z(1), '*k');
% plot3(xE, yE, zE,'LineStyle','--','Color','r','LineWidth',0.3);   % Initial orbit
% plot3(xM, yM, zM,'LineStyle','-.','Color','b','LineWidth',0.3);   % Final orbit
hold on
grid on; 
    
legend('off')
plot3(x, y, z,'k','LineWidth',1);
plot3(x(end), y(end), z(end),'*k');
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1)) * gamma, 'k','LineWidth',1)
plot(tau, u * gamma, 'LineWidth', 0.3)
yline(T * gamma, '--k')
xlabel('Flight time')
ylabel('$\mathbf{a}$')
legend('$u$','$u_x$','$u_y$','$u_z$')
grid on;

figure 
hold on
plot(tau, rad2deg(atan2(u(2,:),u(1,:)))); 
hold off 
grid on;
xlabel('$t$')
ylabel('$\theta$')
title('Thrust in-plane angle')

figure 
hold on
plot(tau, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
hold off 
grid on;
xlabel('$t$')
ylabel('$\phi$')
title('Thrust out-of-plane angle')
