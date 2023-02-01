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

% Initial boundary conditions
initial_state = [0 1];

% Final bounday conditions
final_state = [0 -1]; 

% Spacecraft parameters 
T = 0.5e-4;              % Maximum acceleration 
TOF = 1*365*3600*24;     % Desired TOF for the time-fixed problem

% Setup 
setup.order = n; 
setup.basis = basis;
setup.grid = time_distribution; 
setup.nodes = m; 
setup.FreeTime = true;
setup.resultsFlag = true; 

%% Optimization
% Simple solution    
tic
[C, cost, u, tf, tfapp, tau, exitflag, output] = sb_solver(initial_state, final_state, TOF, L, setup);
toc 

% Compute the optimal analytical solution 
t = (tau+1)/2;
uopt = [-(2/3)/L*(1-t(t <= 3*L)/(3*L)) 0*ones(1,length(t(t > 3*L & t <= 1-3*L))) -(2/3)/L*(1-(1-t(t > 1-3*L & t <= 1))/(3*L))];
    
%% Plot results
% Results 
figure
plot(tau, C); 
legend('$x$', '$\dot{x}$')
xlabel("$t$")
ylabel("$\mathbf{s}$")
grid on;

figure
hold on
plot(tau, uopt);
scatter(tau, u)
xlabel("$t$")
ylabel("$\mathbf{u}$")
grid on;
