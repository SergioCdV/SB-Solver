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
r0 = 6378e3;                      % 1 Earth's radius [m]
mu = 3.986e14;                    % Gavitational parameter of the Earth [m^3 s^−2]
t0 = sqrt(r0^3/mu);               % Fundamental time unit

% Initial orbital elements
coe_initial = [1.5*r0 1e-4]; 
theta0 = deg2rad(110);
coe_initial = [coe_initial theta0]; 

% Final orbital elements
coe_final = [1*r0 0.09]; 
thetaf = deg2rad(260);
coe_final = [coe_final thetaf]; 

%% Initial time of flight
%tfapp = initial_tof(mu, T, initial, final);
tfapp = 100*t0;

%% Normalization
% Gravitational parameter of the body
mu = 1;

% Boundary conditions
coe_initial(1) = coe_initial(1)/r0;
coe_final(1) = coe_final(1)/r0;

initial = [coe_initial zeros(size(coe_initial))];
final = [coe_final zeros(size(coe_final))];

% Sanity check on the SMA evolution 
if (final(1) > initial(1))
    error('Infeasible transfer with Differential Drag');
end

% Time of flight
tfapp = tfapp/t0;

%% Initial approximation to the problem
% Initial guess for the boundary control points
m = 300;    
tau = collocation_grid(m, time_distribution);
[Papp, Capp] = initial_approximation(tau, tfapp, initial, final, basis);

% Initial fitting for n+1 control points
[P0, C0] = initial_fitting(n, tau, Capp, basis);

% Final collocation grid and basis
m = 60;    
tau = collocation_grid(m, time_distribution);
B = state_basis(n, tau, basis);

%% Optimisiation
% Initial guess 
x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
x0 = [x0; tfapp];
L = length(x0)-1;

% Upper and lower bounds (empty in this case)
P_lb = [-Inf*ones(L,1); 0; 0];
P_ub = [Inf*ones(L,1); Inf];

% Objective function
objective = @(x)cost_function(mu, initial, final, n, tau, x, B, basis);

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear constraints
nonlcon = @(x)constraints(mu, initial, final, n, x, B, basis);

% Modification of fmincon optimisation options and parameters (according to the details in the paper)
options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'iter-detailed', 'Algorithm', 'sqp');
options.MaxFunctionEvaluations = 1e6;

% Optimisation
[sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);

% Solution 
P = reshape(sol(1:end-1), [size(P0,1) size(P0,2)]);     % Optimal control points
tf = sol(end);                                          % Optimal time of flight

P(:,[1 2 end-1 end]) = boundary_conditions(tf, n, initial, final, basis);

% Final constraints
[c,ceq] = constraints(mu, initial, final, n, sol, B, basis);

% Final state evolution
C = evaluate_state(P,B,n);
r = sqrt(C(1,:).^2+C(3,:).^2);

% Mass evolution
time = tau*tf;

% Control input
u = acceleration_control(mu,C,tf);
u = u/tf^2;

%% Results
display_results(exitflag, output, t0, tfapp, tf);
plots(); 