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

%% Measurements
% System data 
r0 = 6378e3;                            % 1 Earth's radius [m]
mu = 3.986e14;                          % Gavitational parameter of the Earth [m^3 s^−2]
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

% Time measurements 
t = [0 60 120 180];
 
% Position measurements 
r = zeros(3,length(t)); 
coe_earth = [2*r0 1e-4 0 deg2rad(1) 0 0]; 
theta = [0 pi/4 pi/2 pi];

for i = 1:length(t)
    coe_earth(end) = theta(i);
    s = coe2state(mu, coe_earth);
    r(:,i) = s(1:3);
end 
 
%% Normalization
% Gravitational parameter of the body
mu = 1;

% Measurements normalization 
t = t/t0;       % Time measurements normalization
r = r/r0;       % Position measurements normalization

%% Initial approximation to the problem
% Initial guess for the boundary control points
m = 300;    
tau = collocation_grid(m, time_distribution);
[Papp, Capp] = initial_approximation(tau, t, r, basis); 

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
L = length(x0);

% Upper and lower bounds (empty in this case)
P_lb = -Inf*ones(L,1);
P_ub = Inf*ones(L,1);

% Objective function
objective = @(x)cost_function(mu, initial, final, n, tau, x, B, basis, time_distribution);

% Linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Non-linear constraints
nonlcon = @(x)constraints(mu, T, initial, final, n, x, B, basis, time_distribution);

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
[c,ceq] = constraints(mu, T, initial, final, n, sol, B, basis, time_distribution);

% Final state evolution
C = evaluate_state(P,B,n);
r = sqrt(C(1,:).^2+C(3,:).^2);

% Solution normalization
switch (time_distribution)
    case 'Sundman'
        % Dimensional time grid
        r = sqrt(C(1,:).^2+C(3,:).^2);              % Position vector
        h = angular_momentum(C);                    % Angular momentum vector
        eta = r.^2./h;                              % Mean motion

        % Normalised time grid
        time = tau*tf;

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

%% Results
display_results(exitflag, output, r0, t0, tfapp, tf, dV);
plots(); 