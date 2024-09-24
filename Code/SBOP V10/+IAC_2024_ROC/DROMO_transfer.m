%% Project: SBOPT %%
% Date: 01/08/22

%% 3D low-thrust transfer %% 
% This script provides a main interface to solve 3D low-thrust transfers in DROMO %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 25;                                  % Polynomial order in the state vector expansion
m = 100;                                % Number of sampling points

solver = Solver(basis, n, time_distribution, m);

% Spacecraft parameters 
T = 0.5e-3;              % Maximum acceleration 
N = 0;                   % Number of revolutions

% System data 
r0 = 149597870700;              % 1 AU [m]
mu = 1.32712440042e+20;         % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);             % Fundamental time unit 
Vc = r0/t0;                     % Characteristic velocity

mu = 1;                         % Normalized parameter
gamma = r0/t0^2;                % Characteristic acceleration

%% Boundary conditions 
% Initial orbital elements 
initial_coe = [r0 0 deg2rad(0) deg2rad(0) deg2rad(0) deg2rad(0)];                
theta0 = OrbitalDynamics.KeplerSolver(initial_coe(2), initial_coe(end));
e = initial_coe(2);
sin_theta = sqrt(1-e.^2) .* sin(theta0) ./ (1+e.*cos(theta0));
cos_theta = (cos(theta0)+e) ./ (1+e.*cos(theta0));
E0 = atan2(sin_theta, cos_theta);

if (E0 < 0)
    E0 = E0 + 2*pi;
end

% Final orbital elements 
final_coe = [1.5236*r0 0 deg2rad(0) deg2rad(90) deg2rad(0) deg2rad(270)];  
thetaf = OrbitalDynamics.KeplerSolver(final_coe(2), final_coe(end));
e = final_coe(2);
sin_theta = sqrt(1-e.^2) .* sin(thetaf) ./ (1+e.*cos(thetaf));
cos_theta = (cos(thetaf)-e) ./ (1+e.*cos(thetaf));
Ef = atan2(sin_theta, cos_theta);

if (Ef < 0)
    Ef = Ef + 2*pi;
end

%% Normalization 
initial_coe(1) = initial_coe(1) / r0;   % Normalized initial osculating semi major axis
final_coe(1) = final_coe(1) / r0;       % Normalized final osculating semi major axis
T = T/gamma;                            % Normalized acceleration

%% Boundary conditions (II)
S0 = OrbitalDynamics.coe2dromo(mu, initial_coe);                  % Initial DROMO
SF = OrbitalDynamics.coe2dromo(mu, final_coe);                    % Final DROMO

%% Problem definition 
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 7;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

problem_params = [mu; T; final_coe(2); S0(end); thetaf; N];
S0 = S0(1:7);
SF = SF(1:7);

% Create the problem
OptProblem = IAC_2024_ROC.DROMO(S0, SF, L, StateDimension, ControlDimension, problem_params);

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
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

time = mean(time);

%% Results
% Main plots 
C(4:7,:) = C(4:7,:) ./ sqrt( dot(C(4:7,:), C(4:7,:), 1) );
C = [C(1:7,:); tau; C(8:end,:)];

S = zeros(6,length(tau));
for i = 1:length(tau)
    S(:,i) = OrbitalDynamics.dromo2state(C(1:8,i));
end

% Time of flight
Gamma(1,:) = C(3,:).^3 .* (1 + C(1,:) .* cos(tau) + C(2,:) .* sin(tau)).^2;
deltaT = trapz(tau, 1 ./ Gamma, 2);

%% Plots
% Main plots
x = S(1,:);
y = S(2,:); 
z = S(3,:);
    
% Earth's orbit
thetaE = linspace(0, 2*pi, size(C,2));
s = OrbitalDynamics.coe2state(mu, initial_coe);
initial = OrbitalDynamics.cylindrical2cartesian(s, false).';

% Transfer orbit
s = zeros(6,length(thetaE));

for i = 1:length(thetaE)
    s(:,i) = OrbitalDynamics.coe2state(mu, [initial_coe(1:end-1) initial(2)+thetaE(i)]);
end
xE = s(1,:);
yE = s(2,:);
zE = s(3,:);
    
% Final orbit
s = OrbitalDynamics.coe2state(mu, final_coe);
final = OrbitalDynamics.cylindrical2cartesian(s, false).';

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
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
plot3(0, 0, 0, '*k');
plot3(x(1), y(1), z(1), '*k');
plot3(xE, yE, zE, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 0.3);   
plot3(xM, yM, zM, 'LineStyle', '-.', 'Color', 'b', 'LineWidth', 0.3);   
hold on
grid on; 
legend('off')
plot3(x, y, z, 'k', 'LineWidth', 1);
plot3(x(end), y(end), z(end), '*k');
grid on;
% xticklabels(strrep(xticklabels, '-', '$-$'));
% yticklabels(strrep(yticklabels, '-', '$-$'));
% zticklabels(strrep(zticklabels, '-', '$-$'));

%%

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1)) * gamma, 'k','LineWidth',1)
plot(tau, u * gamma, 'LineWidth', 0.3)
yline(T * gamma, '--k')
xlabel('$s$')
ylabel('$\mathbf{a}$')
legend('$\|\mathbf{a}\|_2$','$a_1$','$a_2$','$a_3$')
grid on;
xlim([min(tau) max(tau)])
yticklabels(strrep(yticklabels, '-', '$-$'));

%%
figure 
hold on
plot(tau, rad2deg(atan2(u(2,:),u(1,:)))); 
hold off 
grid on;
xlabel('$s$')
ylabel('$\theta$')
title('Thrust in-plane angle')

figure 
hold on
plot(tau, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
hold off 
grid on;
xlabel('$s$')
ylabel('$\phi$')
title('Thrust out-of-plane angle')