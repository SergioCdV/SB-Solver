%% Project: 
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif
fig = 1;                                % Figure start number

%% Setup of the collocation method
time_distribution = 'Laguerre';           % Distribution of time intervals
basis = 'Laguerre';                    % Polynomial basis to be use
n = [9 9 9];                            % Order of Bezier curve functions for each coordinate

%% Boundary conditions 
% System data 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

% Earth's orbital elements
coe_earth = [r0 1e-4 0 deg2rad(0) 0]; 
theta0 = deg2rad(95);
coe_earth = [coe_earth theta0]; 

% Mars' orbital elements 
coe_mars = [1.05*r0 1e-4 deg2rad(0) deg2rad(1) deg2rad(0)]; 
thetaf = deg2rad(95);
coe_mars = [coe_mars thetaf]; 

% Initial input revolutions 
K = 0;

% Initial state vector 
s = coe2state(mu, coe_earth);
initial = cylindrical2cartesian(s, false).';

% Final state vector 
s = coe2state(mu, coe_mars);
final = cylindrical2cartesian(s, false).';

%% Initial time of flight
% Spacecraft propulsion parameters 
T = 0.05e-3;     % Maximum acceleration 

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

% Add additional revolutions 
final(2) = final(2)+2*pi*K;

% Time of flight
tfapp = tfapp/t0;

% Spacecraft propulsion parameters 
T = T*(t0^2/r0);

%% Initial approximation to the problem
% Initial guess for the boundary control points
m = 300;   
tau = collocation_grid(m, time_distribution, '');
[Papp, Capp, Napp, tfapp] = initial_approximation(tau, tfapp, initial, final, basis); 

% Initial fitting for n+1 control points
[P0, C0] = initial_fitting(n, tau, Capp, basis);

% Final collocation grid and basis
m = 60;  
tau = collocation_grid(m, time_distribution, '');
[B, tau] = state_basis(n, tau, basis);

%% Optimisiation
% Initial guess 
x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
x0 = [x0; tfapp; Napp];
L = length(x0)-2;

% Upper and lower bounds (empty in this case)
P_lb = [-Inf*ones(L,1); 0; 0];
P_ub = [Inf*ones(L,1); Inf; Inf];

% Objective function
objective = @(x)cost_function(mu, initial, final, n, tau, x, B, basis, time_distribution);

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear constraints
nonlcon = @(x)constraints(mu, T, initial, final, n, x, B, basis, time_distribution, tau);

% Modification of fmincon optimisation options and parameters (according to the details in the paper)
options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'iter-detailed', 'Algorithm', 'sqp');
options.MaxFunctionEvaluations = 1e6;

% Optimisation
[sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);

% Solution 
P = reshape(sol(1:end-2), [size(P0,1) size(P0,2)]);     % Optimal control points
tf = sol(end-1);                                        % Optimal time of flight
N = floor(sol(end));                                    % Optimal number of revolutions 

P = boundary_conditions(tf, n, initial, final, N, P, B, basis);

% Final constraints
[c,ceq] = constraints(mu, T, initial, final, n, sol, B, basis, time_distribution, tau);

% Final state evolution
C = evaluate_state(P,B,n);
r = sqrt(C(1,:).^2+C(3,:).^2);

% Solution normalization
switch (time_distribution)
    case 'Sundman'
        % Normalised time grid
        options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);
        [s,time] = ode45(r,tau,0,options);
        s = tf*s; 
        time = tf*time; 

        % Control input
        u = acceleration_control(mu,C,tf,time_distribution);
        u = u/tf^2;

        % Trajectory cost
        dV = dV/tf;

    otherwise
        % Normalised time grid
        time = tau*tf;

        % Control input
        u = acceleration_control(mu,C,tf,time_distribution);
        u = u/tf^2;

        % Trajectory cost
        dV = dV/tf;
end

% Time domain normalization 
switch (time_distribution)
    case 'Chebyshev'
        tau = (1/2)*(1+tau);
    case 'Legendre'
        tau = (1/2)*(1+tau);
end

%% Results
display_results(exitflag, output, r0, t0, tfapp, tf, dV);
plots(); 
