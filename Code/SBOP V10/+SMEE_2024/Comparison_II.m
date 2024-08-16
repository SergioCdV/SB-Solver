%% Project: SBOPT %%
% Date: 01/08/22
% Date: 29/06/2024

%% Comparison of classical MEE against ideal MEE %% 
% This script provides a comparison in terms of performance between classical Modified Equinoctial Elements and Ideal Stereographic ones %

% The transfer example is provided in Peterson, Arya and Junking, 2023, Connecting the Equinoctial Elements and Rodrigues Parameters: A New Set of
% Elements

%% Set up 
close all
clear

rng(13)

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

% 138925(2001 AU43) orbital elements 
final_mean = [2.83738E11, 0.3765, 2.2567, 1.2593, 2.60614, 0.634857];

% Spacecraft parameters 
T = 5E-4;                       % Maximum acceleration [m/s^2] 
T = T/gamma;                    % Normalized acceleration

% Mission clocks
t0 = 0;                         % Initial normalized clock

L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
ControlDimension = 3;           % Dimension of the control vector

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 10;                                 % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points

solver = Solver(basis, n, time_distribution, m);

%% Optimization
setup.resultsFlag = false; 

% Average results 
iter = 1;                              % Number of revolutions 
time = zeros(3,iter);                   % Convergence time
conv = zeros(3,iter);                   % Convergence flags
feval = zeros(3,iter);                  % Function evaluations
iterations = zeros(3,iter);             % Iterations
cost = zeros(3,iter);                   % Cost function
ToF = zeros(3,iter);                    % Time of Flight
false_positive = zeros(3,iter);         % False positive flag

S0 = zeros(6,iter);
SF = zeros(6,iter);

% Generate the problems 
for i = 1:iter
    % Random boundary conditions 
    initial_coe = initial_mean;
    initial_coe(1) = initial_coe(1) / r0;
    S0(1:6,i) = OrbitalDynamics.coe2equinoctial(initial_coe, true).';       % Initial MEEs

    final_coe = final_mean;
    final_coe(1) = final_coe(1) / r0;
    SF(1:6,i) = OrbitalDynamics.coe2equinoctial(final_coe, true).';         % Final MEEs
end

index = 1; 

% Compute the transfers
if (1)
    % Optimization of the transfers
    for i = index:iter
        % Problem parameters
        problem_params = [mu; T; S0(end,i); SF(end,i); 1 + i];                      
    
        % Regularized motion
        s0 = S0(1:end-1,i);
        sf = SF(1:end-1,i);
    
        % Create the problems 
        StateDimension = 5;                                     % Dimension of the configuration vector. Note the difference with the state vector
    
        OptProblemRMEE = SMEE_2024.LowThrustRMEE(s0, sf, L, StateDimension, ControlDimension, problem_params);
    
        StateDimension = 6;                                     % Dimension of the configuration vector. Note the difference with the state vector
        s0 = [s0; 0];
        sf = [sf; 0];
        OptProblemIMEE = SMEE_2024.LowThrustIMEE(s0, sf, L, StateDimension, ControlDimension, problem_params);

        s0(4:5,1) = s0(4:5,1) ./ ( 1 + sqrt( 1 + dot(s0(4:5,1), s0(4:5,1),1) ) );       % Initial MRP
        sf(4:5,1) = sf(4:5,1) ./ ( 1 + sqrt( 1 + dot(sf(4:5,1), sf(4:5,1),1) ) );       % Final MRP
        OptProblemSIMEE = SMEE_2024.LowThrustSIMEE(s0, sf, L, StateDimension, ControlDimension, problem_params);
    
        % Optimization in classical regularized MEEs
        tic 
        [C_mee, dV, u_mee, Tc, tf, tau_mee, exitflag, output] = solver.solve(OptProblemRMEE);
            
        % Additional check 
        l = tau_mee;
        w = 1 + C_mee(2,:) .* cos(l) + C_mee(3,:) .* sin(l); 
        r = C_mee(1,:) ./ w;
        false_positive(1,i) = ( max(r) > 10 );

        % Results
        time(1,i) = toc;
        conv(1,i) = 1 * (exitflag > 0) + 0 * (exitflag <= 0);
        feval(1,i) = output.funcCount;
        iterations(1,i) = output.iterations;
        cost(1,i) = dV; 
        ToF(1,i) = l(end);

        % Optimization in classical regularized ideal MEEs
        tic 
        [C_imee, dV, u_imee, Tc, tf, tau_imee, exitflag, output] = solver.solve(OptProblemIMEE);
        
        % Additional check 
        l = tau_imee + C_imee(6,:);
        w = 1 + C_imee(2,:) .* cos( l ) + C_imee(3,:) .* sin( l ); 
        r = C_imee(1,:) ./ w;
        false_positive(2,i) = ( max(r) > 10 );

        % Results
        time(2,i) = toc;
        conv(2,i) = 1 * (exitflag > 0) + 0 * (exitflag <= 0);
        feval(2,i) = output.funcCount;
        iterations(2,i) = output.iterations;
        cost(2,i) = dV; 
        ToF(2,i) = tau_imee(1,end) + C_imee(6,end);

        % Optimization in classical regularized stereographic ideal MEEs
%         tic 
%         [C_simee, dV, u_simee, Tc, tf, tau_simee, exitflag, output] = solver.solve(OptProblemSIMEE);
%                 
%         % Additional check 
%         l = tau_simee + C_simee(6,:);
%         w = 1 + C_simee(2,:) .* cos( l ) + C_simee(3,:) .* sin( l ); 
%         r = C_simee(1,:) ./ w;
%         false_positive(3,i) = ( max(r) > 10 );
% 
%         % Results
%         time(3,i) = toc;
%         conv(3,i) = 1 * (exitflag > 0) + 0 * (exitflag <= 0);
%         feval(3,i) = output.funcCount;
%         iterations(3,i) = output.iterations;
%         cost(3,i) = dV; 
%         ToF(3,i) = tau_simee(1,end) + C_simee(6,end);
   
        fprintf('Iteration: %d\n', i);
    end

    % save('+SMEE_2024\rev_comparison.mat')
end

%% Post-processing of results
% Mean values
mu_conv = sum(conv,2) / size(conv,2);
mu_time = mean(time,2);
mu_feval = mean(feval,2);
mu_iter = mean(iterations,2);
mu_cost = mean(cost,2);
mu_tof = mean(ToF,2);

% Compute the different histograms
% GenHistogram(conv(1,:), conv(2,:), conv(3,:), 2, 'r', 'b', 'g', 'Convergence');
% 
% GenHistogram(false_positive(1,:) & conv(1,:), false_positive(2,:) & conv(2,:), false_positive(3,:) & conv(3,:), 2, 'r', 'b', 'g', 'False positive');
% 
% GenHistogram(time(1,:), time(2,:), time(3,:), 20, 'r', 'b', 'g', 'Comp. time [s]');
% GenHistogram(feval(1,:), feval(2,:), feval(3,:), 20, 'r', 'b', 'g', 'Func. evaluations');
% GenHistogram(iterations(1,:), iterations(2,:), iterations(3,:), 20, 'r', 'b', 'g', 'Iterations');
% GenHistogram(cost(1,:), cost(2,:), cost(3,:), 20, 'r', 'b', 'g', '$J$');
% GenHistogram(ToF(1,:), ToF(2,:), ToF(3,:), 20, 'r', 'b', 'g', '$l_f$');

%% Plots
type = 'IMEE';

switch type 
    case 'MEE'
        C = C_mee;
        U = u_mee;
        Tau = tau_mee;
    case 'IMEE' 
        C = C_imee;
        U = u_imee;
        Tau = tau_imee + C_imee(6,:);
    case 'SIMEE'
        C = C_simee;
        U = u_simee;
        Tau = tau_simee + C_simee(6,:);

        % Transformation to RP
        norm_s = dot(C(4:5,:),C(4:5,:),1);
        C(4:5, norm_s > 1) = -C(4:5, norm_s > 1) ./ norm_s(norm_s > 1);  
        norm_s = dot(C(4:5,:),C(4:5,:),1);
        C(4:5,:) = 2 * C(4:5,:) ./ (1 - norm_s);
end

% Main plots 
[S] = OrbitalDynamics.equinoctial2ECI(mu, [C(1:5,:); Tau], true);
x = S(1,:);
y = S(2,:); 
z = S(3,:);
    
% Earth's orbit
thetaE = linspace(0, 2*pi, size(C,2));
    
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
xlabel('$x$ [-]')
ylabel('$y$ [-]')
zlabel('$z$ [-]')
plot3(xE,yE,zE,'LineStyle','--','Color','r','LineWidth',1);   % Earth's orbit
plot3(xM,yM,zM,'LineStyle','-.','Color','b','LineWidth',1);   % Mars' orbit
plot3(x,y,z,'k','LineWidth',0.3);
legend('$\mathcal{C}_0$', '$\mathcal{C}_f$', '$\mathcal{C}(L_{\alpha})$', 'AutoUpdate', 'off')
plot3(0,0,0,'*k');
plot3(x(1),y(1),z(1),'*k');
plot3(x(end),y(end),z(end),'*k');
hold on
grid on; 
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));
% axis('equal')
grid on;

%%

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(Tau, sqrt(dot(U,U,1)) * gamma, 'k','LineWidth',1)
plot(Tau, U * gamma, 'LineWidth', 0.3)
yline(T * gamma, '--k')
xlabel('$l$')
xlim([min(Tau) max(Tau)])
ylabel('$\mathbf{u} [\textrm{m}/\textrm{s}^2]$')
legend('$\|\mathbf{u}\|$', '$u_r$','$u_t$','$u_n$', '$u_{max}$')
grid on;

figure 
hold on
plot(Tau, rad2deg(atan2(U(2,:),U(1,:)))); 
hold off 
grid on;
xlabel('$l$')
xlim([min(Tau) max(Tau)])
ylabel('$\theta$')
title('Thrust in-plane angle')

figure 
hold on
plot(Tau, rad2deg(atan2(U(3,:), sqrt(U(1,:).^2+U(2,:).^2)))); 
hold off 
grid on;
xlabel('$l$')
xlim([min(Tau) max(Tau)])
ylabel('$\phi$')
title('Thrust out-of-plane angle')

% Position coordinates
% figure_coordinates = figure;
% title('Spacecraft position coordinates in time')
% hold on

% subplot(3,1,1)
% hold on
% plot(tau + C(6,:), xM, 'LineStyle','-.','Color','b','LineWidth',0.3)
% plot(tau + C(6,:), xE, 'LineStyle','--','Color','r','LineWidth',0.3)
% plot(tau + C(6,:), x, 'k','LineWidth',1)
% plot(tau(1) + C(6,1), x(1),'*k','DisplayName','')
% plot(tau(end) + C(6,end),x(end),'*k','DisplayName','')
% xlabel('$L_f$')
% xlim([min(tau + C(6,:)) max(tau + C(6,:))])
% ylabel('$x$ [AU]')
% grid on;
% 
% subplot(3,1,2)
% hold on
% plot(tau+ C(6,:), yM, '-.','LineWidth',0.3)
% plot(tau+ C(6,:), yE, '--','LineWidth',0.3)
% plot(tau+ C(6,:), y, 'k','LineWidth',1)
% plot(tau(1) + C(6,1), y(1),'*k','DisplayName','')
% plot(tau(end) + C(6,end),y(end),'*k','DisplayName','')
% xlabel('$L_f$')
% xlim([min(tau + C(6,:)) max(tau + C(6,:))])
% ylabel('$y$ [AU]')
% grid on;
% 
% subplot(3,1,3)
% hold on
% plot(tau + C(6,:), zM, '-.','LineWidth',0.3)
% plot(tau + C(6,:), zE, '--','LineWidth',0.3)
% plot(tau + C(6,:), z, 'k','LineWidth',1)
% plot(tau(1) + C(6,1), z(1),'*k','DisplayName','')
% plot(tau(end) + C(6,end),z(end),'*k','DisplayName','')
% xlabel('$L_f$')
% xlim([min(tau + C(6,:)) max(tau + C(6,:))])
% ylabel('$z$ [AU]')
% grid on;

function GenHistogram(datos1, datos2, datos3, numBins, color1, color2, color3, xlab)
    % generarHistogramaComparativo genera y muestra histogramas comparativos de dos conjuntos de datos
    %
    % Inputs:
    %   datos1   - Vector de datos para el primer histograma
    %   datos2   - Vector de datos para el segundo histograma
    %   datos3   - Vector de datos para el tercer histograma
    %   numBins  - Número de bins (barras) para los histogramas
    %   color1   - Color de las barras del primer histograma (opcional)
    %   color2   - Color de las barras del segundo histograma (opcional)
    %   color3   - Color de las barras del tercer histograma (opcional)
    %
    % Ejemplo de uso:
    %   generarHistogramaComparativo(randn(1000,1), randn(1000,1)*2, 20, 'r', 'b');

    % Validar los inputs
    if nargin < 4
        color1 = 'b'; % Color por defecto para el primer histograma es azul
        color2 = 'r'; % Color por defecto para el segundo histograma es rojo
        color3 = 'g'; % Color por defecto para el tercer histograma es verde
    elseif nargin < 5
        color2 = 'r'; % Color por defecto para el segundo histograma es rojo
    end
    
    % Outliers
    if any( diff(datos1) )
        idx = abs( datos1 - mean(datos1,2) ) < 3 * std(datos1, [], 2);
        datos1 = datos1(1,idx);
    end

    if any( diff(datos2) )
        idx = abs( datos2 - mean(datos2,2) ) < 3 * std(datos2, [], 2);
        datos2 = datos2(1,idx);
    end

    if any( diff(datos3) )
        idx = abs( datos3 - mean(datos3,2) ) < 3 * std(datos3, [], 2);
        datos3 = datos3(1,idx);
    end

    % Crear la figura
    figure;
    hold on; 

    plot(1:size(datos1,2), datos1, 'Color', color1, 'Marker', 'o')
    plot(1:size(datos2,2), datos2, 'Color', color2, 'Marker', 'square')
    plot(1:size(datos3,2), datos3, 'Color', color3, 'Marker', '^')

    ylabel(xlab);
    xlabel('$N$');
    legend('MEE', 'IMEE', 'SIMEE')
    grid on;

    % Personalizar el gráfico
    hold off;
end

