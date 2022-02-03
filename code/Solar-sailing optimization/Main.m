%% Project: 
% Date: 01/02/22

%% Main script
% Version 5 
% Following the paper "Initial design..." by Fan et. al

% This script aims to perform trajectory design optimisation processes
% based on Bernstein polynomials collocation methods

%% Graphicss
animations = 0;     % Set to 1 to generate the gif
fig = 1;            % Figure start number

%% Variables to be defined for each run
m = 200;                                  % Number of discretization points
time_distribution = 'Gauss-Lobatto';     % Distribution of time intervals
sigma = 1;                               % If normal distribution is selected

%% Constraints
amax = 1.5e-4;                           % Maximum acceleration available [m/s^2]

%% Collocation method 
% Order of Bezier curve functions for each coordinate
n = [12 12 16];

%% Global constants
r0 = 149597870700;                      % 1 au [m] (for dimensionalising)
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2]

%% Initial definitions
% Generate the time interval discretization distribution
switch (time_distribution)
    case 'Linear'
        tau = linspace(0,1,m); % non-dimensional time
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
    otherwise
        error('An appropriate time array distribution must be specified')
end

%% Boundary conditions of the problem
% Initial data
[initial, final] = initial_data(r0, 1);

% Initial guess for the boundary control points
[tfapp, Papp, ~, Capp] = initial_approximation(mu, r0, amax, tau, initial, final);

% Initial fitting for n+1 control points
[B, P0, C0] = initial_fitting(n, tau, Capp, 'Non-orthogonal');

%% Optimisiation
% Upper and lower bounds (empty in this case)
P_lb = [];
P_ub = [];

% Objective function
objective = @(P)velocity_variation(mu, r0, tau, tfapp, P, B, n);

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear constraints
nonlcon = @(P)constraints(mu, r0, tfapp, n, P, P0, B, amax);

% Modification of fmincon optimisation options and parameters (according to the details in the paper)
options = optimoptions('fmincon', 'TolFun', 1e-6, 'Display', 'iter-detailed', 'Algorithm', 'sqp');
options.MaxFunctionEvaluations = 5e+03;

% Optimisation
[P, dV, exitflag, output] = fmincon(objective, P0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);

% Dimensionalising 
tf_final = flight_time(P, B, m, tfapp, r0, n);  % Final time of flight
dV = dV*(r0/tfapp);                             % Final velocity change

%% Results
% State vector approximation calculation
C = evaluate_state(P,B,n);

% Results
display_results(P0, P, B, m, exitflag, output, tfapp, r0, n)
plots();
