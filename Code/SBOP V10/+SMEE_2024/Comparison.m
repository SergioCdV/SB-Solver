%% Project: SBOPT %%
% Date: 01/08/22
% Date: 15/05/2024

%% Comparison of classical MEE against ideal MEE %% 
% This script provides a comparison in terms of performance between classical Modified Equinoctial Elements and Ideal Stereographic ones %

% The transfer example is provided in Peterson, Arya and Junking, 2023, Connecting the Equinoctial Elements and Rodrigues Parameters: A New Set of
% Elements

%% Set up 
close all
clear

%% Problem definition 
% System data 
r0 = 149597870700;              % 1 AU [m]
mu = 1.32712440042e+20;         % Gavitational parameter of the Sun [m^3 s^−2] 
Tc = sqrt(r0^3/mu);             % Fundamental time unit 
Vc = r0/Tc;                     % Characteristic velocity
gamma = r0/Tc^2;                % Characteristic acceleration
mu = 1;                         % Normalized parameter

% Earth's orbital elements
initial_mean = [1.497251E11, 0.0173, 2.8152, 7.6438E-05, 5.2940, 0.7221];
initial_sig = [1E3, 0.005, 0.1, 5E-6, 0.1, 0].^2;

% 138925(2001 AU43) orbital elements 
final_mean = [2.83738E11, 0.3765, 2.2567, 1.2593, 2.60614, 0.634857];
final_sig = [1E3, 0.005, 0.1, 0.1, 0.1, 2*pi].^2;

% Spacecraft parameters 
T = 5E-3;                       % Maximum acceleration [m/s^2] 
T = T/gamma;                    % Normalized acceleration

% Mission clocks
t0 = 0;                         % Initial normalized clock

L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
ControlDimension = 3;           % Dimension of the control vector

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 7;                                 % Polynomial order in the state vector expansion
m = 100;                                % Number of sampling points

solver = Solver(basis, n, time_distribution, m);

%% Optimization
% Average results 
iter = 100; 
time = zeros(2,iter);                   % Convergence time
conv = zeros(2,iter);                   % Convergence flags
feval = zeros(2,iter);                  % Function evaluations
iterations = zeros(2,iter);             % Iterations
cost = zeros(2,iter);                   % Cost function
ToF = zeros(2,iter);                    % Time of Flight

setup.resultsFlag = false; 

% Optimization of the transfers
for i = 1:iter
    % Random boundary conditions 
    initial_coe = mvnrnd(initial_mean, diag(initial_sig), 1);
    initial_coe(1) = initial_coe(1) / r0;
    S0 = OrbitalDynamics.coe2equinoctial(initial_coe, true).';       % Initial MEEs

    final_coe = mvnrnd(final_mean, diag(final_sig), 1);
    final_coe(end) = 2*pi * rand;
    final_coe(1) = final_coe(1) / r0;
    SF = OrbitalDynamics.coe2equinoctial(final_coe, true).';         % Final MEEs

    problem_params = [mu; T; S0(end); SF(end)];                      % Problem parameters

    % Regularized motion
    S0 = S0(1:end-1);
    SF = SF(1:end-1);

    % Create the problems 
    StateDimension = 5;                                     % Dimension of the configuration vector. Note the difference with the state vector

    OptProblemRMEE = SMEE_2024.LowThrustRMEE(S0, SF, L, StateDimension, ControlDimension, problem_params);

    StateDimension = 6;                                     % Dimension of the configuration vector. Note the difference with the state vector
    S0 = [S0; 0];
    SF = [SF; 0];
    OptProblemSMEE = SMEE_2024.LowThrustSMEE(S0, SF, L, StateDimension, ControlDimension, problem_params);

    % Optimization in classical MEEs
    tic 
    [C, dV, u, Tc, tf, tau, exitflag, output] = solver.solve(OptProblemRMEE);
    time(1,i) = toc;
    conv(1,i) = 1 * (exitflag > 0) + 0 * (exitflag <= 0);
    feval(1,i) = output.funcCount;
    iterations(1,i) = output.iterations;
    cost(1,i) = dV; 
    ToF(1,i) = tf;

    % Optimization in classical SMEEs
    tic 
    [C, dV, u, Tc, tf, tau, exitflag, output] = solver.solve(OptProblemSMEE);
    time(2,i) = toc;
    conv(2,i) = 1 * (exitflag > 0) + 0 * (exitflag <= 0);
    feval(2,i) = output.funcCount;
    iterations(2,i) = output.iterations;
    cost(2,i) = dV; 
    ToF(2,i) = tf + C(6,end);
end

%% Post-processing of results
% Compute the different histograms
GenHistogram(time(1,:), time(2,:), 100, 'r', 'b', 'Comp. time [s]');
GenHistogram(conv(1,:), conv(2,:), 2, 'r', 'b', 'Convergence');
GenHistogram(feval(1,:), feval(2,:), 100, 'r', 'b', 'Func. evaluations');
GenHistogram(iterations(1,:), iterations(2,:), 100, 'r', 'b', 'Iterations');
GenHistogram(cost(1,:), cost(2,:), 100, 'r', 'b', '$J$');
GenHistogram(ToF(1,:), ToF(2,:), 100, 'r', 'b', '$L_f$');

%% Plots
% Main plots 
[S] = OrbitalDynamics.equinoctial2ECI(mu, [C(1:5,:); C(6,:) + tau], true);
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

function GenHistogram(datos1, datos2, numBins, color1, color2, xlab)
    % generarHistogramaComparativo genera y muestra histogramas comparativos de dos conjuntos de datos
    %
    % Inputs:
    %   datos1   - Vector de datos para el primer histograma
    %   datos2   - Vector de datos para el segundo histograma
    %   numBins  - Número de bins (barras) para los histogramas
    %   color1   - Color de las barras del primer histograma (opcional)
    %   color2   - Color de las barras del segundo histograma (opcional)
    %
    % Ejemplo de uso:
    %   generarHistogramaComparativo(randn(1000,1), randn(1000,1)*2, 20, 'r', 'b');

    % Validar los inputs
    if nargin < 4
        color1 = 'b'; % Color por defecto para el primer histograma es azul
        color2 = 'r'; % Color por defecto para el segundo histograma es rojo
    elseif nargin < 5
        color2 = 'r'; % Color por defecto para el segundo histograma es rojo
    end

    % Crear la figura
    figure;

    % Crear el primer histograma
    idx = abs( datos1 - mean(datos1,2) ) < 3 * std(datos1, [], 2);
    histogram(datos1(1,idx), numBins, 'FaceColor', color1, 'FaceAlpha', 0.5, 'Normalization', 'probability');
    hold on; % Mantener el primer histograma para superponer el segundo

    % Crear el segundo histograma
    idx = abs( datos2 - mean(datos2,2) ) < 3 * std(datos2, [], 2);
    histogram(datos2(1,idx), numBins, 'FaceColor', color2, 'FaceAlpha', 0.5, 'Normalization', 'probability');

    % Personalizar el gráfico
    xlabel(xlab);
    ylabel('Prob.');
    legend('MEE', 'SMEE');
    grid on;
    hold off;
end

