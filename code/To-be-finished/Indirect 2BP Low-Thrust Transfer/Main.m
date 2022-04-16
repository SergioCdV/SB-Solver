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
n = repmat([5 5 5 5 5 5],1,2);

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
% Gravitational parameter of the body
mu = 1; 

% Thruser/accleration and spacecraft mass data
T = 1.405e-1; 
m0 = 1/T; 
Isp = 0.07/T;

% Earth orbital element 
coe_earth = [1 1e-4 0 0 0]; 
s = coe2state(mu, [coe_earth deg2rad(0)]);
initial = cylindrical2cartesian(s, false).';
initial = [initial zeros(size(initial))];

% Mars orbital elements 
coe_mars = [1 1e-4 0 deg2rad(0) 0]; 
s = coe2state(mu, [coe_mars deg2rad(90)]);
final = cylindrical2cartesian(s, false).';
final = [final zeros(size(final))];

% Initial guess for the boundary control points
[Papp, ~, Capp, tfapp] = initial_approximation(mu, tau, n, T, initial, final, 'Bernstein');

% Initial fitting for n+1 control points
[B, P0, C0] = initial_fitting(n, tau, Capp, 'Orthogonal Bernstein');

%% Optimisiation
% Initial guess 
x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
x0 = [x0; tfapp];
L = length(x0)-1;

% Upper and lower bounds (empty in this case)
P_lb = [-Inf*ones(L,1); 0];
P_ub = [Inf*ones(L,1); 10];

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear constraints
nonlcon = @(x)constraints(mu, m0, Isp, T, tau, initial, final, n, m, x, B);

% Objective function
objective = @(x)cost_function(T,x,B,m,n,tau);

% Modification of fmincon optimisation options and parameters (according to the details in the paper)
options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'iter-detailed', 'Algorithm', 'sqp');
options.MaxFunctionEvaluations = 1e6;

% Optimisation
[sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);

% Solution 
[c,ceq] = constraints(mu, m0, Isp, T, tau, initial, final, n, m, sol, B);
P = reshape(sol(1:end-1), [size(P0,1) size(P0,2)]);
tf = sol(end);
C = evaluate_state(P,B,n);
time = tau*tf;

% Dimensionalising
C(13:end,:) = C(13:end,:)/tf;

% Compute the control vector 
u = zeros(3,size(C,2));

%% Results
display_results(exitflag, output, tfapp, tf);
plots(); 