%% Project: SBOPT %%
% Date: 01/08/22
% Date: 29/06/2024

%% Comparison of classical MEE against ideal MEE %% 
% This script provides a comparison in terms of performance between classical Modified Equinoctial Elements and Ideal Stereographic ones %

% The transfer example is given by Conway et al., maximum energy planar transfer

%% Set up 
close all
clear

rng(13)

%% Problem definition 
% System data 
mu = 1;                         % Normalized parameter

% Earth's orbital elements
initial_mean = [1, 0, 0, 0, 0, 0];

% 138925(2001 AU43) orbital elements 
final_mean = [3, 0, 0, 0, 0, 0];

% Spacecraft parameters 
T = 0.01;                       % Normalized acceleration

% Mission clocks
t0 = 0;                         % Initial normalized clock
tf = 50;                        % Final normalized clock
delta_t = tf - t0;              % Time of flight

L = 1;                          % Degree of the dynamics (maximum derivative order of the ODE system)
ControlDimension = 2;           % Dimension of the control vector

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 25;                                  % Polynomial order in the state vector expansion
m = 100;                                % Number of sampling points

solver = Solver(basis, n, time_distribution, m);

%% Optimization
setup.resultsFlag = false; 

% Average results 
iter = 10;                               % Number of revolutions 
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
    S0(1:6,i) = OrbitalDynamics.coe2equinoctial(initial_coe, true).';       % Initial MEEs

    final_coe = final_mean;
    SF(1:6,i) = OrbitalDynamics.coe2equinoctial(final_coe, true).';         % Final MEEs
end 

% Compute the transfers
if (1)
    % Optimization of the transfers
    for i = 1:iter
        % Problem parameters
        problem_params = [mu; T; S0(end,i); SF(end,i); 4; delta_t];                      
    
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
%         tic 
%         [C_mee, dV, u_mee, Tc, tf, tau_mee, exitflag, output] = solver.solve(OptProblemRMEE);
%             
%         % Additional check 
%         l = tau_mee;
%         w = 1 + C_mee(2,:) .* cos(l) + C_mee(3,:) .* sin(l); 
%         r = C_mee(1,:) ./ w;
%         false_positive(1,i) = ( max(r) > 10 );
% 
%         % Results
%         time(1,i) = toc;
%         conv(1,i) = 1 * (exitflag > 0) + 0 * (exitflag <= 0);
%         feval(1,i) = output.funcCount;
%         iterations(1,i) = output.iterations;
%         cost(1,i) = dV; 
%         ToF(1,i) = tf;

        % Optimization in regularized ideal MEEs
%         tic 
%         [C_imee, dV, u_imee, Tc, tf, tau_imee, exitflag, output] = solver.solve(OptProblemIMEE);
%         
%         % Additional check 
%         l = tau_imee + C_imee(6,:);
%         w = 1 + C_imee(2,:) .* cos(l) + C_imee(3,:) .* sin(l); 
%         r = C_imee(1,:) ./ w;
%         false_positive(2,i) = ( max(r) > 10 );
% 
%         % Results
%         time(2,i) = toc;
%         conv(2,i) = 1 * (exitflag > 0) + 0 * (exitflag <= 0);
%         feval(2,i) = output.funcCount;
%         iterations(2,i) = output.iterations;
%         cost(2,i) = dV; 
%         ToF(2,i) = l(end);

        % Optimization in regularized stereographic ideal MEEs
        tic 
        [C_simee, dV, u_simee, Tc, tf, tau_simee, exitflag, output] = solver.solve(OptProblemSIMEE);
                
        % Additional check 
        l = tau_simee + C_simee(6,:);
        w = 1 + C_simee(2,:) .* cos(l) + C_simee(3,:) .* sin(l); 
        r = C_simee(1,:) ./ w;
        false_positive(3,i) = ( max(r) > 10 );

        % Results
        time(3,i) = toc;
        conv(3,i) = 1 * (exitflag > 0) + 0 * (exitflag <= 0);
        feval(3,i) = output.funcCount;
        iterations(3,i) = output.iterations;
        cost(3,i) = dV; 
        ToF(3,i) = l(end);
   
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
% GenHistogram(false_positive(1,:) & conv(1,:), ~false_positive(2,:) & conv(2,:), false_positive(3,:) & conv(3,:), 2, 'r', 'b', 'g', 'False positive');
% 
% GenHistogram(time(1,:), time(2,:), time(3,:), 20, 'r', 'b', 'g', 'Comp. time [s]');
% GenHistogram(feval(1,:), feval(2,:), feval(3,:), 20, 'r', 'b', 'g', 'Func. evaluations');
% GenHistogram(iterations(1,:), iterations(2,:), iterations(3,:), 20, 'r', 'b', 'g', 'Iterations');
% GenHistogram(cost(1,:), cost(2,:), cost(3,:), 20, 'r', 'b', 'g', '$J$');
% GenHistogram(ToF(1,:), ToF(2,:), ToF(3,:), 20, 'r', 'b', 'g', '$l_f$');

%% Plots
type = 'SIMEE';

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
        C(4:5,:) = 2 * C(4:5,:) ./ (1 - dot(C(4:5,:),C(4:5,:),1));
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
final_coe = OrbitalDynamics.coe2equinoctial( [C(1:5,end); Tau(end)], false);
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
legend('$\mathcal{C}_0$', '$\mathcal{C}_f$', '$\mathcal{C}(l)$', 'AutoUpdate', 'off')
plot3(0,0,0,'*k');
plot3(x(1),y(1),z(1),'*k');
plot3(x(end),y(end),z(end),'*k');
hold on
grid on; 
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));
axis("equal")
grid on;

%%
% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(Tau, sqrt(dot(U,U,1)), 'k','LineWidth',1)
plot(Tau, U, 'LineWidth', 0.3)
yline(T, '--k')
xlabel('$l$')
xlim([min(Tau) max(Tau)])
ylabel('$\mathbf{u}$')
legend('$\|\mathbf{u}\|$','$u_r$','$u_t$','$u_n$', '$u_{max}$')
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