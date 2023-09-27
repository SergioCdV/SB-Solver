%% Project: SBOPT %%
% Date: 05/05/23

%% Asteroid Landing %% 
% This script provides a main interface to solve the asteroid landing problem %

%% Set up
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                   % Polynomial basis to be use. Alternatively: Legendre, Bernestein, Orthogonal Bernstein
time_distribution = 'Legendre';       % Distribution of time intervals. Alternatively: Bernstein, Orthogonal Bernstein, Chebsyhev, Legendre, Linear, Newton-Cotes, Normal, Random, Trapezoidal
n = 7;                                % Polynomial order in the state vector expansion
m = 20;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

% Problem parameters 
Lc = 490;                       % Characteristic length
Tc = 3600;                      % Characteristic time
gamma = Lc/Tc^2;                % Characteristic acceleration

T = 20;                         % Maximum control acceleration [m/s^2]
g = 1.6;                        % Maximum acceleration [m/s^2]
a = [1 2 1].';                  % Ellipsoid semiaxis

T = T/Lc;
g = g/Lc;

% Boundary conditions
SF = [0 * a(1); 0 * a(2); a(3); zeros(3,1)];            % Final conditions

problem_params = [g; a; T];

%% Optimization
n_problems = 1;    
k = 1;

% Sampling of the initial conditions
Sigma = diag([0.1 * a(1); 0.1 * a(2); 1 * a(3)]);
r0 = [a(1) * zeros(1,n_problems); a(2) * zeros(1,n_problems); a(3) * ones(1,n_problems)] + ...
     + mvnrnd(zeros(3,1), Sigma, n_problems).';

for i = 1:n_problems
    % Create the problem
    S0 = [r0(:,i); 0; 0; 0];               % Initial conditions
    OptProblem = Problems.AsteroidLanding(S0, SF, L, StateDimension, ControlDimension, problem_params);
    
    tic
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    toc 

    % Save the solution
    if ((exitflag == 1) || (exitflag == 2) || (exitflag == 3))
        State{k} = C; 
        control{k} = u; 
        time{k} = tau;
        k = k+1;
    end
end

% Create the .mat file 
if (exist('State', 'var'))
%     save landing_set State control time T g a;
end

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

%% Plots
% Descend representation
figure;
hold on
plot(tau, C(1:6,:))
xlabel('$t$')
ylabel('$\mathbf{s}$')
legend('$x$', '$y$', '$z$', '$\dot{x}$', '$\dot{y}$', '$\dot{z}$')
hold on
grid on; 
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure_orbits = figure;
hold on
plot(tau, C(4:6,:))
xlabel('$t$')
ylabel('$\dot{\mathbf{r}}$')
hold on
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
yline(T, '--k')
plot(tau, sqrt(dot(u,u,1)), 'k','LineWidth',1)
plot(tau, u, 'LineWidth', 0.3)
xlabel('$t$')
ylabel('$\mathbf{a}$')
legend('$T_{max}$', '$\|\mathbf{u}\|$', '$u_x$', '$u_y$', '$u_z$', 'Location', 'East')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));