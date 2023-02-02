%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Setup of the solution method
% Modify this!!!
time_distribution = 'Legendre';       % Distribution of time intervals (chose between Legendre, Chebyshev or Bernstein)
basis = 'Legendre';                   % Polynomial basis to be use
n = 20;                               % Polynomial order in the state vector expansion
m = 30;                               % Number of sampling points

% Initial boundary conditions
initial_state = [0 0];

% Final bounday conditions
final_state = [1 0]; 

% Constraints 
TOF = 1;                              % Desired TOF for the time-fixed problem

% Setup 
setup.order = n; 
setup.basis = basis;
setup.grid = time_distribution; 
setup.nodes = m; 
setup.FreeTime = false;
setup.resultsFlag = true; 

%% Optimization
% Simple solution    
tic
[C, cost, u, tf, t, exitflag, output] = sb_solver(initial_state, final_state, TOF, setup);
toc 

% Compute the optimal analytical solution (for the original boundary conditions and TOF)
uopt = -12*t+6;
    
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
