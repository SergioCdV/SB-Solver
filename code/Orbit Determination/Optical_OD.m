%% Project: 
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the solution method
time_distribution = 'Linear';     % Distribution of time intervals
basis = 'Chebyshev';              % Polynomial basis to be use
cost_policy = 'Least Squares';    % Minimization cost function
n = [7];                            % Order of Bezier curve functions for each coordinate
m = 100;                          % Number of sampling points

% System data 
Re = 6378e3;                      % Mean Earth radius
r0 = 7000e3;                      % Normalizing orbital distance [m]
mu = 3.986e14;                    % Gavitational parameter of the Earth [m^3 s^−2] 
J2 = 1.08263e-3;                  % Second zonal harmonic of the Earth
t0 = sqrt(r0^3/mu);               % Fundamental time unit

system.mu = mu; 
system.coeffs = J2;
system.distance = r0; 
system.radius = Re;     
system.time = t0; 

% Initial conditions
initial_coe = [r0 1e-3 0 deg2rad(0) 0]; 
theta0 = deg2rad(90);
initial_coe = [initial_coe theta0]; 

% Spacecraft true trajectory 
dt = 1e-3;            % Integration time step
tf = 3600;            % Final integration epoch
tspan = 0:dt:tf;      % Integration time span 

options = odeset('RelTol', 2.25e-14, 'AbsTol', 1e-22);

s0 = coe2state(mu, initial_coe);
u = [0;0;-1e-1];
[t,s] = ode45(@(t,s)dynamics(mu, J2, Re, t, s, u), tspan, s0, options);

% Final conditions 
final_coe = state2coe(mu, s(end,:).', 'Inertial');
final_coe = final_coe(1:end-1);

% Generation of the measurements 
K = [1e3 1e3 1e3 1e3 1e3 1e3];
K = repmat(K, size(s,1), 1);
S = s-K+2*K.*rand(size(s));
measurements = [t S(:,1:3)./sqrt(S(:,1).^2+S(:,2).^2+S(:,3).^2)].'; 

index = randperm(size(measurements,2), 5); 
index = sort(index);
measurements = measurements(:,index);

% Spacecraft propulsion parameters 
T = 1;     % Maximum acceleration 

% Setup 
setup.resultsFlag = true; 
setup.animations = false; 

%% Results
% Simple solution    
tic
[C, e, u, tf, tfapp, tau, exitflag, output] = sbod_optimization(system, tf, initial_coe, final_coe, measurements, T, m, n, cost_policy, time_distribution, basis, setup);
toc 

E = abs(1e-1-abs(u(3,:)));
umean = mean(E); 
umax = max(E); 
usigma = std(E);

figure(1)
hold on
plot3(s(:,1)/r0, s(:,2)/r0, s(:,3)/r0); 
scatter3(S(index,1)/r0, S(index,2)/r0, S(index,3)/r0, 'filled', 'b');

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, e, u, tf, tfapp, tau, exitflag, output] = sbod_optimization(system, tf, initial_coe, final_coe, measurements, T, m, n, cost_policy, time_distribution, basis, setup);
    time(i) = toc;
end

time = mean(time);

