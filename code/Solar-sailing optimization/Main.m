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
m = 90;                                 % Number of discretization points
time_distribution = 'Random';           % Distribution of time intervals
sigma = 1;                              % If normal distribution is selected

%% Constraints
amax = 1.5e-4;                          % Maximum acceleration available [m/s^2]

%% Collocation method 
% Order of Bezier curve functions for each coordinate
n = 20;

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
    otherwise
        error('An appropriate time array distribution must be specified')
end

%% Boundary conditions of the problem
% Initial data
[initial, final] = initial_data(r0, 1);

% Initial guess for the boundary control points
[tfapp, Papp, ~, Capp] = initial_approximation(mu, r0, amax, tau, initial, final);

% Initial fitting for n+1 control points
[B, P0, C0] = initial_fitting(n, tau, Capp);

%% Optimisiation
% Upper and lower bounds (empty in this case)
P_lb = [];
P_ub = [];

% Objective function
objective = @(P)velocity_variation(mu, r0, tau, tfapp, P, B);

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear constraints
nonlcon = @(P)constraints(mu, r0, tfapp, P, P0, B, amax);

% Modification of fmincon optimisation options and parameters (according to the details in the paper)
options = optimoptions('fmincon','TolCon',1e-9,'TolFun',1e-6,'Display','iter-detailed');
options.MaxFunctionEvaluations = 5e+03;

% Optimisation
[P, dV, exitflag, output] = fmincon(objective, P0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);

% Dimensionalising 
tf_final = flight_time(P, B, m, tfapp, r0);     % Final time of flight
dV = dV*(r0/tfapp);                             % Final velocity change

%% Results
% Coordinate calculation
C = zeros(9,size(B,2));     % Preallocation for speed
k = size(B,1)/3;            % Number of control points
for i = 1:3
    C(1+3*(i-1):3*i,:) = P*B(1+k*(i-1):k*i,:);
end

display_results(P0, P, B, m, exitflag, output, tfapp, r0, n)
plots();
