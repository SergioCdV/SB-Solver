%% Project: SBOPT %%
% Date: 05/05/23

%% Kalman LQR %% 
% This script provides a main interface to solve Kalman's LQR problem %

%% Set up
close all
clear
rng(1)

%% Numerical solver definition 
basis = 'Legendre';                   % Polynomial basis to be use. Alternatively: Legendre, Bernestein, Orthogonal Bernstein
time_distribution = 'Legendre';       % Distribution of time intervals. Alternatively: Bernstein, Orthogonal Bernstein, Chebsyhev, Legendre, Linear, Newton-Cotes, Normal, Random, Trapezoidal
n = [7 7 7];                          % Polynomial order in the state vector expansion
m = 100;                              % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

% Boundary conditions
SF = zeros(6,1);                % Final conditions

% Problem parameters
tf = 20; 
A = rand(3,6);
problem_params = [tf; reshape(A, [], 1)];

%% Optimization
n_problems = 1;    
k = 1;

% Sampling of the initial conditions
Sigma = eye(6);
r0 = mvnrnd(zeros(6,1), Sigma, n_problems).';
TOF = 20 + 80 * rand(1,n_problems);
%%
tic
for i = 1:n_problems
    % Create the problem
    problem_params(1) = TOF(i);
    S0 = r0(:,i);                  % Initial conditions
    OptProblem = Problems.LQR(S0, SF, L, StateDimension, ControlDimension, problem_params);
    
    tic
    [C, dV, u, t0, tf, tau, exitflag, output, P] = solver.solve(OptProblem);
    toc 

    % Save the solution
    if ((exitflag == 1) || (exitflag == 2) || (exitflag == 3))
        State{k} = C; 
        control{k} = u; 
        time{k} = tau;
        k = k+1;
    end
end
toc

% Create the .mat file 
% if (exist('State', 'var'))
%     save LQR_set State control time;
% end

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output, P] = solver.solve(OptProblem);
    time(i) = toc;
end

time = mean(time);
%%
M1 = PolynomialBases.Bezier().LB_tmatrix(7);
M2 = PolynomialBases.Bezier().DB_tmatrix(7, 1);

B =  M1 * PolynomialBases.Bezier().basis(7, tau/tf);
dB = (M1 * M2.') * PolynomialBases.Bezier().basis(7, tau/tf);

L = PolynomialBases.Legendre().derivative(7, 2 * tau/tf-1, 1);

C2(1:3,:) = P * B; 
C2(4:6,:) = P * dB; 

%% Plots
figure;
hold on
plot(tau, C(1:3,:))
plot(tau, C2(1:3,:), '--k')
xlabel('$t$')
ylabel('$\mathbf{s}$')
hold on
grid on; 
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));

figure_orbits = figure;
hold on
plot(tau, C(4:6,:))
plot(tau, C2(4:6,:), '--k')
xlabel('Flight time')
ylabel('$\dot{\mathbf{r}}$')
hold on
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1)), 'k','LineWidth',1)
plot(tau, u, 'LineWidth', 0.3)
xlabel('$t$')
ylabel('$\mathbf{a}$')
legend('$\|u\|$', '$u_1$', '$u_2$', '$u_3$')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
