%% Project: 
% Date: 01/02/22

%% Main script
% Version 5 
% Following the paper "Initial design..." by Fan et. al

% This script aims to perform trajectory design optimisation processes
% based on Bernstein polynomials collocation methods

%% Graphics
set_graphics(); 

animations = 0;     % Set to 1 to generate the gif
fig = 1;            % Figure start number

%% Variables to be defined for each run
m = 50;                                 % Number of discretization points
time_distribution = 'Linear';           % Distribution of time intervals
sigma = 1;                              % If normal distribution is selected

%% Collocation method 
% Order of Bezier curve functions for each coordinate
n = [5 5 5 5];

%% Initial definitions
% Generate the time interval discretization distribution
switch (time_distribution)
    case 'Linear'
        tau = linspace(0,1,m);
    case 'Normal'
        pd = makedist('Normal');
        pd.sigma = sigma;
        xpd = linspace(-3,3,m);
        tau = cdf(pd,xpd);
    case 'Random'
        tau = rand(1, m);
        tau = sort(tau);
    case 'Gauss-Lobatto'
        i = 1:m;
        tau = -cos((i-1)/(m-1)*pi);
        tau = (tau-tau(1))/(tau(end)-tau(1));
    case 'Legendre-Gauss'
        tau = LG_nodes(0,1,m);
    case 'Bezier'
        tau = B_nodes(0,1,m);
    case 'Orthonormal Bezier'
        tau = OB_nodes(0,1,m);
    otherwise
        error('An appropriate time array distribution must be specified')
end

%% Boundary conditions of the problem
% Initial data
mu = 1; 
r0 = 1; 
rf = 2;
tf = 30;
T = 0.1405; 
m0 = 1/T; 
Isp = 0.5328825;
initial = [r0 0 0 sqrt(mu/r0)]; 
final = [rf pi 0 sqrt(mu/rf)];

% Initial guess for the boundary control points
[Papp, ~, Capp] = initial_approximation(mu, tf, tau, n, initial, final, 'Bernstein');

% Initial fitting for n+1 control points
[B, P0, C0] = initial_fitting(n, tau, Capp, 'Orthogonal Bernstein');

%% Optimisiation
% Initial guess 
x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
x0 = [x0; zeros(2*m,1); tf];
L = length(x0)-1;

% Upper and lower bounds (empty in this case)
P_lb = [-Inf*ones(L,1); 0];
P_ub = [Inf*ones(L,1); Inf];

% Objective function
objective = @(x)cost_function(x,B,m,n,tau);

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear constraints
nonlcon = @(x)constraints(mu, m0, Isp, T, tau, initial, final, n, m, x, B);

% Modification of fmincon optimisation options and parameters (according to the details in the paper)
options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'iter-detailed', 'Algorithm', 'sqp');
options.MaxFunctionEvaluations = 1e6;

% Optimisation
[sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);

% Solution 
[c,ceq] = constraints(mu, m0, Isp, T, tau, initial, final, n, m, sol, B);
P = reshape(sol(1:end-1-2*m), [size(P0,1) size(P0,2)]);
tf = sol(end);
C = evaluate_state(P,B,n);
u = reshape(sol(end-2*m:end-1), [2 m]);
time = tau*tf;

% Dimensionalising
C(4:6,:) = C(4:6,:)/tf;

%% Results
figure 
hold on
plot(time, C(1:4,:)); 
hold off 
grid on;
legend('$r$', '$u$', '$v$')
xlabel('Time')
ylabel('State')
title('State evolution in time')

figure 
hold on
plot(time, u); 
hold off 
grid on;
legend('$u_x$', '$u_y$')
xlabel('Time')
ylabel('$\mathbf{u}$')
title('Acceleration vector')

figure 
hold on
plot(time, unwrap(atan2(u(2,:),u(1,:)))); 
hold off 
grid on;
xlabel('Time')
ylabel('$\theta$')
title('Thrust angle')

figure 
R(1,:) = rf*cos(0:1e-2:2*pi);
R(2,:) = rf*sin(0:1e-2:2*pi);
R(3,:) = r0*cos(0:1e-2:2*pi);
R(4,:) = r0*sin(0:1e-2:2*pi);
r = C(1,:).*[cos(C(2,:)); sin(C(2,:))];
hold on
plot(R(1,:),R(2,:), 'b.-');
plot(r(1,:),r(2,:))
plot(R(3,:),R(4,:), 'r.-');
hold off
legend('Initial orbit', 'Final orbit')
title('Low thrust transfer');
xlabel('$x$ coordinate')
ylabel('$y$ coordinate')
grid on;