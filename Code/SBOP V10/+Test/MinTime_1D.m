%% Project: SBOPT %%
% Date: 01/08/22

%% Minimum time 1D %% 
% This script provides a main interface to solve the 1D minimum time problem %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                   % Polynomial basis to be use
time_distribution = 'Legendre';       % Distribution of time intervals
n = 30;                               % Polynomial order in the state vector expansion
m = 100;                              % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 1;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 1;           % Dimension of the control vector

% Add boundary conditions
S0 = [0 0].';                   % Initial conditions
SF = [1 0].';                   % Final conditions

% Problem parameters 
T = 1;                          % Maximum acceleration

% Create the problem
OptProblem = Problems.MinTime_1D(S0, SF, L, StateDimension, ControlDimension, T);

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

%% Analytical solution 
% Analytical solution
s = (S0-SF)/T;           % Reduced state

% In-plane motion
F = s(1)+sign(s(2))*s(2)^2/2;
Sigma = sign(F);

if (Sigma == 0)
    tf(2) = abs(s(1));
    tspan = linspace(0,tf(2),length(tau));
    uopt = -sign(s)*ones(1,length(tspan));
    A = [s(1)+s(2).*tspan+uopt(1,:).*tspan.^2/2; s(2)+uopt.*tspan];
else
    Lambda = sqrt(Sigma*s(1)+s(2)^2/2);
    Delta = [Lambda+Sigma*s(2); Lambda];
    tf(2) = sum(Delta);
    tspan = linspace(0,tf(2),length(tau));

    % Preallocation 
    A = zeros(2,length(tspan));

    % State trajectory and control
    U = -Sigma*ones(1,length(tspan(tspan < Delta(1))));
    A(:,1:length(tspan(tspan < Delta(1)))) = [s(1)+s(2).*tspan(tspan < Delta(1))+U.*tspan(tspan < Delta(1)).^2/2; s(2)+U.*tspan(tspan < Delta(1))];

    U = Sigma*ones(1,length(tspan(tspan >= Delta(1))));
    Ss = [0.5*(s(1)+0.5*(Sigma)*s(2)^2); -Sigma*Lambda];
    A(:,length(tspan(tspan < Delta(1)))+1:length(tspan)) = [Ss(1)+Ss(2).*(tspan(tspan >= Delta(1))-Delta(1))+U.*(tspan(tspan >= Delta(1))-Delta(1)).^2/2; ...
                                                            Ss(2)+U.*(tspan(tspan >= Delta(1))-Delta(1))];

    % Complete control law 
    uopt = [-Sigma*ones(1,length(tspan(tspan < Delta(1)))) Sigma*ones(1,length(tspan(tspan >= Delta(1))))];
end

%% Plots
% State representation
figure_orbits = figure;
hold on
xlabel('$t$')
ylabel('$\mathbf{s}$')
plot(tau, C(1:2,:));
hold off
legend('off')
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, u, 'LineWidth', 0.3)
plot(tau, uopt, 'LineWidth', 0.3)
xlabel('$t$')
ylabel('$\mathbf{u}$')
grid on;

% Phase space plot 
figure_propulsion = figure;
hold on
plot(C(1,:), C(2,:), 'LineWidth', 0.3)
xlabel('$\mathbf{x}$')
ylabel('$\dot{\mathbf{x}}$')
grid on;