%% Project: 
% Date: 01/02/22

%% Set up
set_graphics(); 
close all

animations = 0;                         % Set to 1 to generate the gif

%% Setup of the solution method
time_distribution = 'Regularized';      % Distribution of time intervals
basis = 'Bernstein';                    % Polynomial basis to be use
n = 10;                                 % Order of Bezier curve functions for each coordinate
m = 60;                                 % Number of sampling points

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

% Earth's orbital elements
initial_coe = [r0 0 deg2rad(0) deg2rad(0) deg2rad(0)]; 
theta0 = deg2rad(0);
initial_coe = [initial_coe theta0]; 
Q0 = [cos(theta0) sin(theta0) 0; -sin(theta0) cos(theta0) 0; 0 0 1];

%% Monte-carlo analysis 
cases = 1e3;                            % Number of Monte Carlo cases to consider

% Orbital elements' statistical distribution 
% a = 0.95+0.1*rand(1,cases);            % Semimajor axis distribution 
% e = 1e-4+1e-3*rand(1,cases);           % Eccentricity distribution
% RAAN = deg2rad(15*rand(1,cases));      % RAAN distribution
% I = deg2rad(5*rand(1,cases));          % Inclination distribution
% omega = deg2rad(15*rand(1,cases));     % AoP distribution
% thetaf = deg2rad(360*rand(1,cases));   % True anomaly distribution 

% Preallocation for speed
% Def_factor = zeros(2, cases);
% V = zeros(cases, cases);
% Tf = zeros(cases, cases);
% S = zeros(cases, cases);

d = zeros(1,cases);
alpha = zeros(1,cases);
psi = zeros(1,cases);
V = zeros(1,cases);
Tf = zeros(1,cases); 
S = zeros(1,cases);

% Orbit statistical distribution 
a = 0.95+0.1*rand(1,cases);            % Semimajor axis distribution 
e = 1e-4+1e-3*rand(1,cases);           % Eccentricity distribution
RAAN = deg2rad(15*rand(1,cases));      % RAAN distribution
I = deg2rad(5*rand(1,cases));          % Inclination distribution
omega = deg2rad(15*rand(1,cases));     % AoP distribution
thetaf = deg2rad(360*rand(1,cases));   % True anomaly distribution 

for i = 1:cases 
    % Final orbital elements
    final_coe = [a(i)*r0 e(i) RAAN(i) I(i) omega(i) thetaf(i)];

    % Orbit difference
    da = initial_coe(1)-a(i)*r0;
    db = initial_coe(1)*sqrt(1-initial_coe(2)^2)-a(i)*r0*sqrt(1-e(i)^2);
    Qf = euler_matrix(final_coe);                          
    alpha(i) = norm([da db]) / initial_coe(1);
    psi(i) = acos((1/2)*(trace(Q0*Qf.')-1));
    d(i) = abs(alpha(i)+1i*psi(i));
   
    % Optimisation
    [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, time_distribution, basis, n, setup);

    % Optimization success
    S(i) = (exitflag == 1) || (exitflag == 2) ||(exitflag == 3);

    % Results
    V(i) = dV*S(i);                    % Final dV cost
    Tf(i) = tf*S(i);                   % Final time of flight
end

% [Def_factor(1,:), index] = sort(Def_factor(1,:)); 
% a = a(index); 
% e = e(index);

% tic 
% for i = 1:cases
%     % Final COE 
%     final_coe = [a(i)*r0 e(i)];
% 
%     for j = 1:cases
%         % Final orbital elements
%         final_coe(3:6) = [RAAN(j) I(j) omega(j) thetaf(j)];
%         Qf = euler_matrix(final_coe);                          % Final Euler matrix
%         Def_factor(2,j) = acos((1/2)*(trace(Q0*Qf.')-1));
% 
%         % Optimisation
%         [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, time_distribution, basis, n, setup);
% 
%         % Optimization success
%         S(i,j) = (exitflag == 1) || (exitflag == 2) ||(exitflag == 3);
% 
%         % Results
%         V(i,j) = dV;                    % Final dV cost
%         Tf(i,j) = tf;                   % Final time of flight
%     end
% end
% toc
% 
% [Def_factor(2,:), index] = sort(Def_factor(2,:)); 
% S = S(:,index);
% V = V(:,index);
% Tf = Tf(:,index);

%% Results
[d, index] = sort(d); 
V = V(index);
Tf = Tf(index); 
S = S(index);

% Convergence plot
figure 
hold on
plot(d, S, 'o')
hold off
grid on; 
ylabel('Convergence $\sigma$')
xlabel('Orbit global relative difference $d$')

figure 
subplot(1,2,1)
plot(d, V, 'o')
grid on; 
ylabel('Normalized transfer cost $\Delta V$')
xlabel('Orbit global relative difference $d$')

subplot(1,2,2)
plot(d, Tf, 'o')
grid on; 
ylabel('Normalized time of flight $t_f$')
xlabel('Orbit global relative difference $d$')


[alpha, index] = sort(alpha); 
V = V(index);
Tf = Tf(index); 
S = S(index);
% Convergence plot
figure 
hold on
plot(alpha, S, 'o')
hold off
grid on; 
ylabel('Convergence $\sigma$')
xlabel('Orbit global relative difference $d$')

figure 
subplot(1,2,1)
plot(alpha, V, 'o')
grid on; 
ylabel('Normalized transfer cost $\Delta V$')
xlabel('Orbit shape difference $\alpha$')

subplot(1,2,2)
plot(alpha, Tf, 'o')
grid on; 
ylabel('Normalized time of flight $t_f$')
xlabel('Orbit shape difference $\alpha$')

[psi, index] = sort(psi); 
V = V(index);
Tf = Tf(index); 
S = S(index);

% Convergence plot
figure 
hold on
plot(psi, S, 'o')
hold off
grid on; 
ylabel('Convergence $\sigma$')
xlabel('Orbit orientation difference $\psi$')

figure 
subplot(1,2,1)
plot(psi, V, 'o')
grid on; 
ylabel('Transfer cost $\Delta V$')
xlabel('Orbit orientation difference $\psi$')

subplot(1,2,2)
plot(psi, Tf, 'o')
grid on; 
ylabel('Time of flight $t_f$')
xlabel('Orbit orientation difference $\psi$')

Emin = min(d); 
Emax = max(d); 
E = linspace(Emin, Emax, 25); 
f = zeros(length(E),1); 
r = zeros(length(E),1);
for i = 1:length(S)
    for j = 1:length(E)-1
        if (d(i) < E(j+1))
            f(j) = f(j)+S(i);
            r(j) = r(j)+1;
            break; 
        end
    end
end

f = f./r;

figure 
bar(E,f*100);
xlabel('Orbit global relative difference $d$')
ylabel('Ratio of converged transfers')
grid on; 

% % Porckchop plot 
% Meshgrid 
% [alpha, psi] = meshgrid(Def_factor(1,:), Def_factor(2,:)); 
% rho = 1e3; 
% labels = num2cell(string(Def_factor(2,:)));

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
