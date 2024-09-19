%% Project: SBOPT %%
% Date: 01/08/22

%% KS low-thrust transfer %% 
% This script provides a main interface to solve 3D low-thrust transfers in KS coordinates %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 10;                                  % Polynomial order in the state vector expansion
m = 100;                                % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 4;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 4;           % Dimension of the control vector

% System data 
r0 = 149597870700;              % 1 AU [m]
mu = 1.32712440042e+20;         % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);             % Fundamental time unit 
Vc = r0/t0;                     % Characteristic velocity

mu = 1;                         % Normalized parameter
gamma = r0/t0^2;                % Characteristic acceleration

% Boundary conditions
initial_coe = [r0 1e-3 0 deg2rad(0) deg2rad(0)];                % Initial orbital elements
theta0 = deg2rad(0);
initial_coe = [initial_coe theta0]; 
initial_coe(1) = initial_coe(1) / r0;
e = initial_coe(2);
sin_theta = sqrt(1-e.^2) .* sin(theta0) ./ (1+e.*cos(theta0));
cos_theta = (cos(theta0)+e) ./ (1+e.*cos(theta0));
E0 = atan2(sin_theta, cos_theta);

% Final orbital elements 
final_coe = [2.2*r0 0.542 deg2rad(82.2) deg2rad(13.6) deg2rad(204.2)];   
thetaf = deg2rad(114.4232);

final_coe = [final_coe thetaf]; 
final_coe(1) = final_coe(1) / r0;
e = final_coe(2);
sin_theta = sqrt(1-e.^2) .* sin(thetaf) ./ (1+e.*cos(thetaf));
cos_theta = (cos(thetaf)-e) ./ (1+e.*cos(thetaf));
Ef = atan2(sin_theta, cos_theta);

S0 = OrbitalDynamics.coe2state(mu, initial_coe);
SF =  OrbitalDynamics.coe2state(mu, final_coe);

S0 = LegoKS.state_mapping(S0(1:6), true, "Ecc"); 
SF = LegoKS.state_mapping(SF(1:6), true, "Ecc");

% Spacecraft parameters 
T = 0.5e-3;              % Maximum acceleration 
T = T/gamma;             % Normalized acceleration

problem_params = [mu; T; E0; Ef; 5];

% Definition of the problem
OptProblem = IAC_2024_ROC.EccKSTransfer(S0, SF, L, StateDimension, ControlDimension, problem_params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Dimensions
u = u ./ dot(C(1:4,:), C(1:4,:), 1).^2; 
u = u(1:3,:);

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

time = mean(time);

dV = dV * Vc;

%% Plots
% Compute the true trajectory in Cartesian space 
S = LegoKS.state_mapping(C, false, "1");

% Main plots 
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
    s(:,i) = OrbitalDynamics.coe2state(mu, [initial_coe(1:end-1) thetaE(i)]);
end
xE = s(1,:);
yE = s(2,:);
zE = s(3,:);
    
% Mars's orbit
for i = 1:length(thetaE)
    s(:,i) = OrbitalDynamics.coe2state(mu, [final_coe(1:end-1) thetaE(i)]);
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
plot3(xE,yE,zE,'LineStyle','--','Color','r','LineWidth',0.3);   
plot3(xM,yM,zM,'LineStyle','-.','Color','b','LineWidth',0.3);   
hold on
grid on; 
    
legend('off')
plot3(x,y,z,'k','LineWidth',1);
plot3(x(end),y(end),z(end),'*k');
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1)) * gamma, 'k','LineWidth',1)
plot(tau, u * gamma, 'LineWidth', 0.3)
yline(T * gamma, '--k')
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