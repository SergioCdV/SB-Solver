%% Project: SBOPT %%
% Date: 01/08/22
% Date: 26/04/2024

%% 3D low-thrust transfer %% 
% This script provides a main interface to solve 3D low-thrust transfers in classical regularized MEE coordinates %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 25;                                  % Polynomial order in the state vector expansion
m = 1500;                                % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);
        
% Spacecraft parameters 
T = 0.5E-3;                     % Maximum acceleration 

% System data 
r0 = 149597870700;              % 1 AU [m]
mu = 1.32712440042e+20;         % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);             % Fundamental time unit 
Vc = r0/t0;                     % Characteristic velocity

mu = 1;                         % Normalized parameter
gamma = r0/t0^2;                % Characteristic acceleration

%% Boundary conditions
% Initial orbital elements
initial_coe = [r0 0.015 deg2rad(183.121) deg2rad(0.004) deg2rad(281.94) deg2rad(0)];                
theta0 = OrbitalDynamics.KeplerSolver(initial_coe(2), initial_coe(end));
e = initial_coe(2);
sin_theta = sqrt(1-e.^2) .* sin(theta0) ./ (1+e.*cos(theta0));
cos_theta = (cos(theta0)+e) ./ (1+e.*cos(theta0));
E0 = atan2(sin_theta, cos_theta);

if (E0 < 0)
    E0 = E0 + 2*pi;
end

% Final orbital elements 
final_coe = [1.5236*r0 0.0934 deg2rad(49.4723) deg2rad(1.8464) deg2rad(286.7860)];
thetaf = deg2rad(355.2065);
e = final_coe(2);
sin_theta = sqrt(1-e.^2) .* sin(thetaf) ./ (1+e.*cos(thetaf));
cos_theta = (cos(thetaf)-e) ./ (1+e.*cos(thetaf));
Ef = atan2(sin_theta, cos_theta);
Mf = Ef - e * sin(Ef);

final_coe = [final_coe Mf]; 

% Final orbital elements 
final_coe = [2.2*r0 0.542 deg2rad(82.2) deg2rad(13.6) deg2rad(204.2) deg2rad(114.4232)];  
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
S0(1:6,1) = OrbitalDynamics.coe2equinoctial(initial_coe, true).';     % Initial MEEs
SF(1:6,1) = OrbitalDynamics.coe2equinoctial(final_coe, true).';       % Final MEEs

%% Problem definition
L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 5;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

% Problem parameters
N = 35;
problem_params = [mu; T; S0(end,1); SF(end,1); N];                      

% Regularized motion
S0 = S0(1:end-1,1);
SF = SF(1:end-1,1);

% Create the problem
OptProblem = SMEE_2024.LowThrustRMEE(S0, SF, L, StateDimension, ControlDimension, problem_params);

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

%% Results
% Compute the true trajectory in Cartesian space 
C = [C(1:5,:); tau; C(6:end,:)];
S = OrbitalDynamics.equinoctial2ECI(mu, C, true);

Gamma(1,:) = dot(C(1:4,:), C(1:4,:), 1);

% Time of flight
deltaT(1) = trapz(tau, Gamma(1,:), 2);

%% Plots
% Main plots 
x = S(1,:);
y = S(2,:); 
z = S(3,:);
    
% Initial orbit
thetaE = linspace(0, 2*pi, size(C,2));
s = OrbitalDynamics.coe2state(mu, initial_coe);
initial = OrbitalDynamics.cylindrical2cartesian(s, false).';

% Transfer orbit
s = zeros(6,length(thetaE));

for i = 1:length(thetaE)
    s(:,i) = OrbitalDynamics.coe2state(mu, [initial_coe(1:end-1) thetaE(i)]);
end

xE = s(1,:);
yE = s(2,:);
zE = s(3,:);
    
% Final orbit
s = OrbitalDynamics.coe2state(mu, final_coe);
final = OrbitalDynamics.cylindrical2cartesian(s, false).';

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
plot3(x, y, z, 'k--', 'LineWidth', 0.4);
plot3(x(end), y(end), z(end), '*k');
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

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

%% GIF animation
N = size(S, 2); 

figure;
view([45 45])
hold on;
grid on;
% axis equal;
filename = 'KSDyonisus.gif';

% xlim([-5 5])
% ylim([-5 5])
% zlim([-5 5])

plot3(0, 0, 0, '*k');
plot3(x(1), y(1), z(1), '*k');
plot3(xE, yE, zE, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 0.3);   
plot3(xM, yM, zM, 'LineStyle', '-.', 'Color', 'b', 'LineWidth', 0.3);   
hold on
grid on; 
legend('off')
plot3(x(end), y(end), z(end), '*k');

xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));


trajectory1 = plot3(x, y, z, 'k', 'LineWidth', 0.2);
point1 = plot3(S(1,1), S(2,1), S(3,1), 'ko', 'MarkerFaceColor', 'k');

% Crear la animación
for t = 1:N
    % Actualizar las trayectorias (pasado)
    set(trajectory1, 'XData', S(1,1:t), 'YData', S(2,1:t), 'ZData', S(3,1:t));
    
    % Actualizar las posiciones actuales (puntos)
    set(point1, 'XData', S(1,t), 'YData', S(2,t), 'ZData', S(3,t));
    
    % Dibujar el cuadro
    drawnow;
    
    % Capturar el fotograma para el GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Escribir en el archivo GIF
    if t == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end