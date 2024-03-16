%% Project: SBOPT %%
% Date: 09/03/24

%% 3D low-thrust transfer in the co-orbital CR3BP problem %% 
% This script provides a main interface to solve 3D low-thrust transfers in the co-orbital CR3BP %

%% Set up 
close all
clear

%% Input data
load COOMOT_2024.mat



%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 15;                                % Polynomial order in the state vector expansion
m = 600;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
Ld = 2;                         % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 2;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

% Initial and final clocks 
t0 = 0;                         % Initial clock 
tf = 2 * pi;                    % Final clock
    
% Initial and final conditions 
rt = target_orbit.Trajectory;
ct = chaser_orbit.Trajectory;

S0 = ct(1,1:6)-rt(1,1:6);                               % Relative initial particle

% Initial spiral conditions 
c2 = cn(2);
alpha(1) = (c2-2+sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
alpha(2) = (c2-2-sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
lambda = sqrt(alpha(1));                                % Hyperbolic unstable eigenvalue
omega(1) = sqrt(-alpha(2));                             % Hyperbolic stable eigenvalue
omega(2) = sqrt(c2);                                    % Center manifold eigenvalue
k = (omega(1)^2+1+2*c2)/(2*omega(1));                   % Contraint on the planar amplitude   

psi = pi/2;                               % Out-of-plane angle
phi = atan2(S0(2), -k * S0(1));           % In-plane angle

A = [ abs(S0(1) / cos(phi)); S0(3) / sin(psi)];

S0 = [A; -(S0(4)-A(1)*omega(1)*sin(phi)) / cos(phi); (S0(6)-A(2)*omega(2)*cos(psi)) / sin(psi)];
SF = zeros(Ld*StateDimension,1);                         % Final particle
    
% Mission parameters 
T = 0.5e-3;                     % Maximum acceleration 
T = T/Gamma;                    % Normalized acceleration

problem_params = [mu; t0; tf; omega.'; k; psi; phi; T];

% Create the problem
OptProblem = Problems.IntegralsTransferCR3BP(S0, SF, Ld, StateDimension, ControlDimension, problem_params);

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

% Average results 
iter = 0; 
time = zeros(1,iter);
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
    time(i) = toc;
end

time = mean(time);

dV = dV * Vc;

%% Results 
% Configuration space coordinates
x =    -C(1,:) .* cos(omega(2)*tau + phi);    % Synodic x coordinate
y = k * C(1,:) .* sin(omega(2)*tau + phi);    % Synodic y coordinate
z =     C(2,:) .* sin(omega(1)*tau + psi);    % Synodic z coordinate

C(1,:) = x; 
C(2,:) = y; 
C(3,:) = z; 

%% Plots
% Orbit representation
figure_orbits = figure;
view(3)
hold on
xlabel('$X$')
ylabel('$Y$')
zlabel('$Z$')
plot3(0,0,0, '*k');
plot3(x(1), y(1), z(1), '*k');
% plot3(xE, yE, zE,'LineStyle','--','Color','r','LineWidth',0.3);   % Initial orbit
% plot3(xM, yM, zM,'LineStyle','-.','Color','b','LineWidth',0.3);   % Final orbit
hold on
grid on; 
    
legend('off')
plot3(x, y, z,'k','LineWidth',1);
plot3(x(end), y(end), z(end),'*k');
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, sqrt(dot(u,u,1)) * gamma, 'k','LineWidth',1)
plot(tau, u * gamma, 'LineWidth', 0.3)
yline(T * gamma, '--k')
xlabel('Flight time')
ylabel('$\mathbf{a}$')
legend('$u$','$u_x$','$u_y$','$u_z$')
grid on;

figure 
hold on
plot(tau, rad2deg(atan2(u(2,:),u(1,:)))); 
hold off 
grid on;
xlabel('$t$')
ylabel('$\theta$')
title('Thrust in-plane angle')

figure 
hold on
plot(tau, rad2deg(atan2(u(3,:),sqrt(u(1,:).^2+u(2,:).^2)))); 
hold off 
grid on;
xlabel('$t$')
ylabel('$\phi$')
title('Thrust out-of-plane angle')

%% Save results in a common file 
save COOMOT_2024_results.mat