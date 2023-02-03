%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
setup_path();
set_graphics(); 
close all

%% Setup of the solution method
time_distribution = 'Legendre';        % Distribution of time intervals (Legendre or Chebyshev)
m = 100;                               % Number of sampling points

% Initial boundary conditions
initial_state = [0 1];

% Final bounday conditions
final_state = [0 -1]; 

% Constraints 
L = 1/6;                 % Upper bound on the x dimension
TOF = 1;                 % Desired TOF for the time-fixed problem

% Setup 
setup.grid = time_distribution; 
setup.nodes = m; 

%% Optimization
% Simple solution    
tic
[C, cost, u, tf, t, exitflag, output] = sb_solver(initial_state, final_state, TOF, L, setup);
toc 

% Compute the optimal analytical solution 
uopt = [-(2/3)/L*(1-t(t <= 3*L)/(3*L)) 0*ones(1,length(t(t > 3*L & t <= 1-3*L))) -(2/3)/L*(1-(1-t(t > 1-3*L & t <= 1))/(3*L))];
    
%% Plot results
% Results 
figure
plot(t, C(1:2,:)); 
legend('$x$', '$\dot{x}$')
xlabel("$t$")
ylabel("$\mathbf{s}$")
grid on;

figure
hold on
plot(t, uopt);
scatter(t, u)
xlabel("$t$")
ylabel("$\mathbf{u}$")
legend('Optimal analytical solution', 'Optimal SBOPT solution')
grid on;
