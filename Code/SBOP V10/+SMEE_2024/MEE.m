%% Project: SBOPT %%
% Date: 01/08/22
% Date: 26/04/2024

%% 3D low-thrust transfer %% 
% This script provides a main interface to solve 3D low-thrust transfers in classical MEE coordinates, no regularization %

% The transfer example is provided in Peterson, Arya and Junking, 2023, Connecting the Equinoctial Elements and Rodrigues Parameters: A New Set of
% Elements

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 7;                                 % Polynomial order in the state vector expansion
m = 200;                                % Number of sampling points

solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 7;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

% System data 
r0 = 149597870700;              % 1 AU [m]
mu = 1.32712440042e+20;         % Gavitational parameter of the Sun [m^3 s^−2] 
Tc = sqrt(r0^3/mu);             % Fundamental time unit 
Vc = r0/Tc;                     % Characteristic velocity

mu = 1;                         % Normalized parameter
gamma = r0/Tc^2;                % Characteristic acceleration

% Earth's orbital elements
initial_coe = [r0 1e-3 0 deg2rad(0) deg2rad(0)]; 
theta0 = deg2rad(0);
initial_coe = [initial_coe theta0]; 

initial_coe = [1.497251E11, 0.0173, 2.8152, 7.6438E-05, 5.2940, 0.7221];

initial_coe(1) = initial_coe(1) / r0;
S0 = OrbitalDynamics.coe2equinoctial(initial_coe, true).';       % Initial MEEs

% Initial mass 
m0 = 2800;          % Initial mass [kg]
S0 = [S0; m0];

% Mars' orbital elements 
final_coe = [1.1*r0 1e-3 deg2rad(0) deg2rad(60) deg2rad(0)]; 
thetaf = deg2rad(100);
final_coe = [final_coe thetaf];

final_coe = [2.83738E11, 0.3765, 2.2567, 1.2593, 2.60614, 0.634857];

final_coe(1) = final_coe(1) / r0;
SF = OrbitalDynamics.coe2equinoctial(final_coe, true).';         % Final MEEs

% Final mass 
m0 = 2800;          % Final mass [kg]
SF = [SF; m0];

% Spacecraft parameters 
g0 = 9.81;               % Reference gravity acceleration [m / s^2]
Isp = 3000;              % Spacecraft specific impulse [s]

Isp = Isp / Tc;          % Normalized Isp
g0 = g0 / gamma;         % Normalized g0

T = 0.45;                % Maximum thrust [N] 
T = T/gamma;             % Normalized acceleration

c = Isp * g0;            % Exhaust velocity

% Mission clocks
t0 = 0;                  % Initial normalized clock
tf = 1720 * 86400 / Tc;  % Final normalized clock

% Problem parameters
problem_params = [mu; T; t0; tf; c; final_coe(end)];

% Create the problem
OptProblem = SMEE_2024.LTransferMEE(S0, SF, L, StateDimension, ControlDimension, problem_params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, Tc, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

time = mean(time);

%% Plots
% Main plots 
[S] = OrbitalDynamics.equinoctial2ECI(mu, C(1:6,:), true);
x = S(1,:);
y = S(2,:); 
z = S(3,:);
    
% Earth's orbit
thetaE = linspace(0, 2*pi, size(C,2));

s = OrbitalDynamics.coe2state(mu, initial_coe);
initial = OrbitalDynamics.cylindrical2cartesian(s, false).';

s = OrbitalDynamics.coe2state(mu, final_coe);
final = OrbitalDynamics.cylindrical2cartesian(s, false).';
    
s = zeros(6,length(thetaE));
for i = 1:length(thetaE)
    s(:,i) = OrbitalDynamics.coe2state(mu, [initial_coe(1:end-1) initial(2)+thetaE(i)]);
end
xE = s(1,:);
yE = s(2,:);
zE = s(3,:);
    
% Mars's orbit
for i = 1:length(thetaE)
    s(:,i) = OrbitalDynamics.coe2state(mu, [final_coe(1:end-1) final(2)+thetaE(i)]);
end
xM = s(1,:);
yM = s(2,:);
zM = s(3,:);

% Orbit representation
figure_orbits = figure;
view(3)
hold on
xlabel('$X$ coordinate')
ylabel('$Y$ coordinate')
zlabel('$Z$ coordinate')
plot3(0,0,0,'*k');
plot3(x(1),y(1),z(1),'*k');
plot3(xE,yE,zE,'LineStyle','--','Color','r','LineWidth',0.3);   % Earth's orbit
plot3(xM,yM,zM,'LineStyle','-.','Color','b','LineWidth',0.3);   % Mars' orbit
hold on
grid on; 
    
legend('off')
plot3(x,y,z,'k','LineWidth',1);
plot3(x(end),y(end),z(end),'*k');
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1))*gamma, 'k','LineWidth',1)
plot(tau, u, 'LineWidth', 0.3)
yline(T, '--k')
xlabel('Flight time')
ylabel('$\mathbf{a}$')
legend('$a$','$a_\rho$','$a_\theta$','$a_z$')
grid on;

figure 
hold on
plot(tau, rad2deg(atan2(u(2,:),u(1,:)))); 
hold off 
grid on;
xlabel('Time')
ylabel('$\theta$')
title('Thrust in-plane angle')

figure 
hold on
plot(tau, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
hold off 
grid on;
xlabel('Time')
ylabel('$\phi$')
title('Thrust out-of-plane angle')

% Position coordinates
figure_coordinates = figure;
title('Spacecraft position coordinates in time')
hold on

subplot(3,1,1)
hold on
plot(tau, xM, 'LineStyle','-.','Color','b','LineWidth',0.3)
plot(tau, xE, 'LineStyle','--','Color','r','LineWidth',0.3)
plot(tau, x, 'k','LineWidth',1)
plot(tau(1), x(1),'*k','DisplayName','')
plot(tau(end),x(end),'*k','DisplayName','')
xlabel('Flight time [days]')
ylabel('$x$ [AU]')
grid on;

subplot(3,1,2)
hold on
plot(tau, yM, '-.','LineWidth',0.3)
plot(tau, yE, '--','LineWidth',0.3)
plot(tau, y, 'k','LineWidth',1)
plot(tau(1), y(1),'*k','DisplayName','')
plot(tau(end),y(end),'*k','DisplayName','')
xlabel('Flight time [days]')
ylabel('$y$ [AU]')
grid on;

subplot(3,1,3)
hold on
plot(tau, zM, '-.','LineWidth',0.3)
plot(tau, zE, '--','LineWidth',0.3)
plot(tau, z, 'k','LineWidth',1)
plot(tau(1), z(1),'*k','DisplayName','')
plot(tau(end),z(end),'*k','DisplayName','')
xlabel('Flight time [days]')
ylabel('$z$ [AU]')
grid on;