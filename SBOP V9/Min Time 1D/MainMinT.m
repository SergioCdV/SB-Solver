%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Problem definition 
% Numerical solver definition 
time_distribution = 'Legendre';       % Distribution of time intervals
basis = 'Legendre';                   % Polynomial basis to be use
n = 30;                                % Polynomial order in the state vector expansion
m = 300;                                % Number of sampling points
T = 1;                                 % Degree of the dynamics 

OptProblem = Problem().DefineSolver(n, basis, m, time_distribution).AddDynamics(length(n), 1, 2); 

% Add boundary conditions
S0 = [0 0].';
SF = [1 0].';
OptProblem = OptProblem.AddBoundaryConditions(S0, SF).AddParameters(T);

% Add functions 
OptProblem = OptProblem.AddFunctions(@(initial, final, beta, t0, tf)BoundaryConditions_minT(initial, final, beta, t0, tf), @(params, beta, t0, tf, tau, s)ControlFunction_minT(params, beta, t0, tf, tau, s), ...
                                     @(params, beta, t0, tf, s, u)CostFunction_minT(params, beta, t0, tf, s, u), @(beta, P)LinConstraints_minT(beta, P), ...
                                     @(params, beta, t0, tf, tau, s, u)NlinConstraints_minT(params, beta, t0, tf, tau, s, u), ...
                                     @BoundsFunction_minT, ...
                                     @(params, initial, final)InitialGuess_minT(params, initial, final));

%% Optimization
% Simple solution    
tic
[C, cost, u, t0, tf(1), tau, exitflag, output] = sb_solver(OptProblem);
toc 

%% Analytical solution 
% Analytical solution
s = (S0-SF)/T;           % Reduced state

% In-plane motion
F = s(1)+sign(s(2))*s(2)^2/2;
Sigma = sign(F);

if (Sigma == 0)
    tf(2) = abs(s(1));
    tspan = linspace(0,tf(2),length(tau));
    uopt = -sign(s)*ones(1,length(tspan));
    A = [s(1)+s(2).*tspan+uopt(1,:).*tspan.^2/2; s(2)+uopt.*tspan];
else
    Lambda = sqrt(Sigma*s(1)+s(2)^2/2);
    Delta = [Lambda+Sigma*s(2); Lambda];
    tf(2) = sum(Delta);
    tspan = linspace(0,tf(2),length(tau));

    % Preallocation 
    A = zeros(2,length(tspan));

    % State trajectory and control
    U = -Sigma*ones(1,length(tspan(tspan < Delta(1))));
    A(:,1:length(tspan(tspan < Delta(1)))) = [s(1)+s(2).*tspan(tspan < Delta(1))+U.*tspan(tspan < Delta(1)).^2/2; s(2)+U.*tspan(tspan < Delta(1))];

    U = Sigma*ones(1,length(tspan(tspan >= Delta(1))));
    Ss = [0.5*(s(1)+0.5*(Sigma)*s(2)^2); -Sigma*Lambda];
    A(:,length(tspan(tspan < Delta(1)))+1:length(tspan)) = [Ss(1)+Ss(2).*(tspan(tspan >= Delta(1))-Delta(1))+U.*(tspan(tspan >= Delta(1))-Delta(1)).^2/2; ...
                                                            Ss(2)+U.*(tspan(tspan >= Delta(1))-Delta(1))];

    % Complete control law 
    uopt = [-Sigma*ones(1,length(tspan(tspan < Delta(1)))) Sigma*ones(1,length(tspan(tspan >= Delta(1))))];
end

%% Plots
% State representation
figure_orbits = figure;
hold on
xlabel('$t$')
ylabel('$\mathbf{s}$')
plot(tau, C(1:2,:));
hold off
legend('off')
grid on;

% Propulsive acceleration plot
figure_propulsion = figure;
hold on
plot(tau, u, 'LineWidth', 0.3)
plot(tau, uopt, 'ok');
xlabel('$t$')
ylabel('$\mathbf{u}$')
grid on;

% Phase space plot 
figure_propulsion = figure;
hold on
plot(C(1,:), C(2,:), 'LineWidth', 0.3)
xlabel('$\mathbf{x}$')
ylabel('$\dot{\mathbf{x}}$')
grid on;