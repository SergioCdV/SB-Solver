%% Project: 
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif
fig = 1;                                % Figure start number

%% Setup of the collocation method
time_distribution = 'Linear';           % Distribution of time intervals
basis = 'Bernstein';                    % Polynomial basis to be use
sigma = 1;                              % If normal distribution is selected
n = [9 9 9];                            % Order of Bezier curve functions for each coordinate

%% Boundary conditions 
% System data 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2]
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

% Earth's orbital elements
coe_earth = [r0 1e-4 0 deg2rad(0) 0]; 
theta0 = deg2rad(110);
coe_earth = [coe_earth theta0]; 

% Mars' orbital elements 
coe_mars = [1.5*r0 0.09 deg2rad(0) deg2rad(0) 0]; 
thetaf = deg2rad(260);
coe_mars = [coe_mars thetaf]; 

% Initial state vector 
s = coe2state(mu, coe_earth);
initial = cylindrical2cartesian(s, false).';

% Final state vector 
s = coe2state(mu, coe_mars);
final = cylindrical2cartesian(s, false).';

%% Initial time of flight
% Spacecraft propulsion parameters 
T = 0.5e-3;     % Maximum acceleration 

% Initial TOF
tfapp = initial_tof(mu, T, initial, final);

%% Normalization
% Gravitational parameter of the body
mu = 1;

% Boundary conditions
coe_earth(1) = coe_earth(1)/r0;
coe_mars(1) = coe_mars(1)/r0;

s = coe2state(mu, coe_earth);
initial = cylindrical2cartesian(s, false).';

s = coe2state(mu, coe_mars);
final = cylindrical2cartesian(s, false).';

% Time of flight
tfapp = tfapp/t0;

% Spacecraft propulsion parameters 
T = T*(t0^2/r0);

%% Initial approximation to the problem
% Initial guess for the boundary control points
m = 300;    
tau = collocation_grid(m, time_distribution);
[Papp, Capp, Napp] = initial_approximation(tau, tfapp, initial, final, basis);

% New initial TOF
tfapp = tfapp*Napp;

% Initial fitting for n+1 control points
%basis = 'Orthogonal Bernstein';
[P0, C0] = initial_fitting(n, tau, Capp, basis);

% Final collocation grid and basis
m = 60;    
tau = collocation_grid(m, time_distribution);
B = state_basis(n, tau, basis);

%% Optimisiation
% Initial guess 
x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
x0 = [x0; tfapp; Napp];
L = length(x0)-2;

% Upper and lower bounds (empty in this case)
P_lb = [-Inf*ones(L,1); 0; 0];
P_ub = [Inf*ones(L,1); Inf; Inf];

% Objective function
objective = @(x)cost_function(mu, initial, final, n, tau, x, B, basis);

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear constraints
nonlcon = @(x)constraints(mu, T, initial, final, n, x, B, basis);

% Modification of fmincon optimisation options and parameters (according to the details in the paper)
options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'iter-detailed', 'Algorithm', 'sqp');
options.MaxFunctionEvaluations = 1e6;

% Optimisation
[sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);

% Solution 
P = reshape(sol(1:end-2), [size(P0,1) size(P0,2)]);     % Optimal control points
tf = sol(end-1);                                        % Optimal time of flight
N = floor(sol(end));                                    % Optimal number of revolutions 

P(:,[1 2 end-1 end]) = boundary_conditions(tf, n, initial, final, N, basis);

% Final constraints
[c,ceq] = constraints(mu, T, initial, final, n, sol, B, basis);

% Final state evolution
C = evaluate_state(P,B,n);
r = sqrt(C(1,:).^2+C(3,:).^2);

% Mass evolution
time = tau*tf;

% Control input
u = acceleration_control(mu,C,tf);
u = u/tf^2;

%% Results
display_results(exitflag, output, r0, t0, tfapp, tf, dV);
plots(); 