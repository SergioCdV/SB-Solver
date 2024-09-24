%% Project: SBOPT %%
% Date: 21/09/24

%% Comparison KS %% 
% This script provides a main interface to compare the results obtained from the KS and ecc. KS formulations %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Bernstein';                    % Polynomial basis to be use
time_distribution = 'Bernstein';        % Distribution of time intervals
n = 7;                                  % Polynomial order in the state vector expansion
m = 100;                                % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);
        
% Spacecraft parameters 
T = 0.5E-3;                     % Maximum acceleration 
N = 1;                          % Number of revolutions

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
final_coe = [1.5236*r0 0 deg2rad(0) deg2rad(0) deg2rad(0) deg2rad(270)];  
thetaf = OrbitalDynamics.KeplerSolver(final_coe(2), final_coe(end));
e = final_coe(2);
sin_theta = sqrt(1-e.^2) .* sin(thetaf) ./ (1+e.*cos(thetaf));
cos_theta = (cos(thetaf)-e) ./ (1+e.*cos(thetaf));
Ef = atan2(sin_theta, cos_theta);

% Final orbital elements 
% final_coe = [1.5236*r0 0.0934 deg2rad(49.4723) deg2rad(1.8464) deg2rad(286.7860)];
% thetaf = deg2rad(355.2065);
% e = final_coe(2);
% sin_theta = sqrt(1-e.^2) .* sin(thetaf) ./ (1+e.*cos(thetaf));
% cos_theta = (cos(thetaf)-e) ./ (1+e.*cos(thetaf));
% Ef = atan2(sin_theta, cos_theta);
% Mf = Ef - e * sin(Ef);
% 
% final_coe = [final_coe Mf]; 

if (Ef < 0)
    Ef = Ef + 2*pi;
end

%% Normalization 
initial_coe(1) = initial_coe(1) / r0;   % Normalized initial osculating semi major axis
final_coe(1) = final_coe(1) / r0;       % Normalized final osculating semi major axis
T = T/gamma;                            % Normalized acceleration

%% Boundary conditions (II)
S0cart = OrbitalDynamics.coe2state(mu, initial_coe);
SFcart = OrbitalDynamics.coe2state(mu, final_coe);

%% Problem definition
% Classical KS
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 4;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 4;           % Dimension of the control vector

S0 = LegoKS.KS_mapping(S0cart(1:6), true, "1", mu); 
SF = LegoKS.KS_mapping(SFcart(1:6), true, "1", mu);

problem_params = [mu; T; 2*pi * N];
    
% Definition of the problem
OptProblem = IAC_2024_ROC.KSTransfer(S0, SF, L, StateDimension, ControlDimension, problem_params);

% Eccentric KS
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 4;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 4;           % Dimension of the control vector

S0 = LegoKS.KS_mapping(S0cart(1:6), true, "Ecc", mu); 
SF = LegoKS.KS_mapping(SFcart(1:6), true, "Ecc", mu);

problem_params = [mu; T; E0; Ef; N];

% Definition of the problem
OptProblem_ecc = IAC_2024_ROC.EccKSTransfer(S0, SF, L, StateDimension, ControlDimension, problem_params);
  
%% Optimization
% Simple solution    
tic
[C, dV_cl, u, t0, tf, tau, exitflag, output] = solver.solve(OptProblem);
toc 

dV(1) = dV_cl;
%%
tic
[C_ecc, dV_ecc, u_ecc, t0_ecc, tf_ecc, tau_ecc, exitflag, output] = solver.solve(OptProblem_ecc);
toc 

dV(2) = dV_ecc;

%% Results
% Classical KS
u = u ./ dot(C(1:4,:), C(1:4,:), 1).^2;     % Control vector
% u = u(1:3,:);

% Compute the true trajectory in Cartesian space 
S = LegoKS.KS_mapping(C, false, "1", mu);

Gamma(1,:) = dot(C(1:4,:), C(1:4,:), 1);

% Eccentric KS
u_ecc = u_ecc ./ dot(C_ecc(1:4,:), C_ecc(1:4,:), 1).^2;     % Control vector
% u_ecc = u_ecc(1:3,:);

% Compute the true trajectory in Cartesian space 
S_ecc = LegoKS.KS_mapping(C_ecc, false, "Ecc", mu);

[~, h] = LegoKS.OscEnergy(mu, C_ecc, "Ecc");
Gamma(2,:) = dot(C_ecc(1:4,:), C_ecc(1:4,:), 1) ./ sqrt(mu .* h);

% Time of flight
deltaT(1) = trapz(tau, Gamma(1,:), 2);
deltaT(2) = trapz(tau_ecc, Gamma(2,:), 2);

% Bilinear error
l(1,:) = LegoKS.bilinear_function( C(1:4,:), C(5:8,:) );
l(2,:) = LegoKS.bilinear_function( C_ecc(1:4,:), C_ecc(5:8,:) );

%% Error 
error = S - S_ecc; 

S = S_ecc;

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