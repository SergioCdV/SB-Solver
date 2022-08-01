%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the solution method
% System data 
r0 = 149597870700;                      % 1 AU [m]
mu = 1.32712440042e+20;                 % Gavitational parameter of the Sun [m^3 s^−2] 
t0 = sqrt(r0^3/mu);                     % Fundamental time unit

system.mu = mu; 
system.distance = r0; 
system.time = t0; 

% Spacecraft propulsion parameters 
T = 0.05e-3;     % Maximum acceleration 

% Initial input revolutions 
K = 1;

% Setup 
setup.resultsFlag = false; 
setup.animations = false; 

% Initial orbital elements
initial_coe = [r0 0 deg2rad(0) deg2rad(0) deg2rad(0)]; 
theta0 = deg2rad(0);
initial_coe = [initial_coe theta0]; 
Q0 = [cos(theta0) sin(theta0) 0; -sin(theta0) cos(theta0) 0; 0 0 1];

% Final orbital elements 
final_coe = [1.05*r0 5e-3 deg2rad(15) deg2rad(1) deg2rad(10)]; 
thetaf = deg2rad(270);
final_coe = [final_coe thetaf]; 

%% Algorithm construction analysis 
% Set up
n = 4:10;                               % Order of Bezier curve functions for each coordinate
m = 30:5:60;                            % Number of sampling points

% Polynomial basis to be use
basis = {'Bernstein', 'Orthogonal Bernstein', 'Chebyshev', 'Legendre'};   

% Sampling grid to be used 
time_distribution = {'Linear', 'Random', 'Normal', 'Regularized', 'Chebyshev'};      

% Preallocation
DV = zeros(length(basis), length(time_distribution), length(n), length(m));     % Final velocity cost
TF = DV;                                                                        % Final time of flight
O = cell(size(DV));                                                             % Optimization output
E = DV;                                                                         % Exitflag

% Analysis
for i = 1:length(basis) 
    for j = length(time_distribution)
        for k = 1:length(n)
            for l = 1:length(m)
                % Optimisation
                [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m(l), char(time_distribution(j)), char(basis(i)), n(k), setup);
            
                % Results
                E(i,j,k,l) = (exitflag == 1) || (exitflag == 2) ||(exitflag == 3);
                DV(i,j,k,l) = dV;           % Final dV cost
                TF(i,j,k,l) = tf;            % Final time of flight
                O{i,j,k,l} = output;        % Optimization output
            end
        end
    end
end

%% Results
% [d, index] = sort(d); 
% V = V(index);
% Tf = Tf(index); 
% S = S(index);
% 
% Convergence plot
% figure 
% hold on
% plot(d, S, 'o')
% hold off
% grid on; 
% ylabel('Convergence $\sigma$')
% xlabel('Orbit global relative difference $d$')
% 
% figure 
% subplot(1,2,1)
% plot(d, V, 'o')
% grid on; 
% ylabel('Normalized transfer cost $\Delta V$')
% xlabel('Orbit global relative difference $d$')
% 
% subplot(1,2,2)
% plot(d, Tf, 'o')
% grid on; 
% ylabel('Normalized time of flight $t_f$')
% xlabel('Orbit global relative difference $d$')
% 
% 
% [alpha, index] = sort(alpha); 
% V = V(index);
% Tf = Tf(index); 
% S = S(index);
% Convergence plot
% figure 
% hold on
% plot(alpha, S, 'o')
% hold off
% grid on; 
% ylabel('Convergence $\sigma$')
% xlabel('Orbit global relative difference $d$')
% 
% figure 
% subplot(1,2,1)
% plot(alpha, V, 'o')
% grid on; 
% ylabel('Normalized transfer cost $\Delta V$')
% xlabel('Orbit shape difference $\alpha$')
% 
% subplot(1,2,2)
% plot(alpha, Tf, 'o')
% grid on; 
% ylabel('Normalized time of flight $t_f$')
% xlabel('Orbit shape difference $\alpha$')
% 
% [psi, index] = sort(psi); 
% V = V(index);
% Tf = Tf(index); 
% S = S(index);
% 
% Convergence plot
% figure 
% hold on
% plot(psi, S, 'o')
% hold off
% grid on; 
% ylabel('Convergence $\sigma$')
% xlabel('Orbit orientation difference $\psi$')
% 
% figure 
% subplot(1,2,1)
% plot(psi, V, 'o')
% grid on; 
% ylabel('Transfer cost $\Delta V$')
% xlabel('Orbit orientation difference $\psi$')
% 
% subplot(1,2,2)
% plot(psi, Tf, 'o')
% grid on; 
% ylabel('Time of flight $t_f$')
% xlabel('Orbit orientation difference $\psi$')
% 
% Emin = min(d); 
% Emax = max(d); 
% E = linspace(Emin, Emax, 25); 
% f = zeros(length(E),1); 
% r = zeros(length(E),1);
% for i = 1:length(S)
%     for j = 1:length(E)-1
%         if (d(i) < E(j+1))
%             f(j) = f(j)+S(i);
%             r(j) = r(j)+1;
%             break; 
%         end
%     end
% end
% 
% f = f./r;
% 
% figure 
% bar(E,f*100);
% xlabel('Orbit global relative difference $d$')
% ylabel('Ratio of converged transfers')
% grid on; 
% 
% % Porckchop plot 
% Meshgrid 
% [alpha, psi] = meshgrid(Def_factor(1,:), Def_factor(2,:)); 
% rho = 1e3; 
% labels = num2cell(string(Def_factor(2,:)));
% 
% figure 
% subplot(1,2,1)
% contour(alpha, psi, V, rho)
% grid on; 
% xlabel('Orbit size difference $\alpha$')
% ylabel('Orbit orientation difference $\psi$')
% 
% subplot(1,2,2)
% contour(alpha, psi, Tf, rho)
% grid on; 
% xlabel('Orbit size difference $\alpha$')
% ylabel('Orbit orientation difference $\psi$')
% colorbar; 
