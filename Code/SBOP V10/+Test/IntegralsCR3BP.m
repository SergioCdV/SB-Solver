%% Project: SBOPT %%
% Date: 09/03/24

%% 3D low-thrust transfer in the CR3BP problem %% 
% This script provides a main interface to solve 3D low-thrust transfers in the CR3BP through integrals of motion %

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
n = 15;                                % Polynomial order in the state vector expansion
m = 200;                               % Number of sampling points
 
solver = Solver(basis, n, time_distribution, m);

%% Problem definition 
L = 2;                          % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 2;             % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;           % Dimension of the control vector

% System data 
mu = 1.2150581623433623E-2;     % Gravitational parameter of the Earth
Lem = 384400e3;                 % Mean distance from the Earth to the Moon
T0 = 28 * 86400;                % Characteristic time of the Earth-Moon system
n = 2*pi / T0;                  % Characteristic frequency of the system
Vc = Lem * n;                   % Characteristic velocity of the Earth-Moon system
gamma = Lem * n^2;              % Characteristic acceleration of the Earth-Moon system

idx = 1;                                                % Libration point selector
Lc = libration_points(mu);                              % Location of the libration points
cn = legendre_coefficients(mu, idx, Lc(end,idx), 2);    % Compute the Legendre coefficients
c2 = cn(2);                                             % Linear dynamics
alpha(1) = (c2-2+sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
alpha(2) = (c2-2-sqrt(9*c2^2-8*c2))/2;                  % Characteristic root
lambda = sqrt(alpha(1));                                % Hyperbolic unstable eigenvalue
omega(1) = sqrt(-alpha(2));                             % Hyperbolic stable eigenvalue
omega(2) = sqrt(c2);                                    % Center manifold eigenvalue
k = (omega(1)^2+1+2*c2)/(2*omega(1));                   % Contraint on the planar amplitude                          

% Initial and final clocks 
t0 = 0;                         % Initial clock 
tf = 2 * pi;                    % Final clock
    
% Initial and final conditions 
S0 = [0.9261 0 0.3616 0 -0.0544  0].';    % State vector of a vertical orbit                 
SF = [1.0406 0 0.1735 0 -0.0770 0].';     % State vector of a butterfly orbit

psi = pi/2;                               % Out-of-plane angle
phi = atan2(S0(2), -k * S0(1));           % In-plane angle

S0 = [5E-3; 1E-3; 0; 0];                  % Initial conditions on the integrals of motion
SF = zeros(4,1);                          % Final conditions on the integrals of motion
    
% Mission parameters 
T = 0.5e-3;                     % Maximum acceleration 
T = T/gamma;                    % Normalized acceleration

problem_params = [mu; t0; tf; omega.'; k; psi; phi; T];

% Create the problem
OptProblem = Problems.IntegralsTransferCR3BP(S0, SF, L, StateDimension, ControlDimension, problem_params);

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

dV = dV * Vc;

%% Results 
% Configuration space coordinates
I = C(1:2,:);         % Integrals of motion
    
%% Plots
% Orbit representation
% figure_orbits = figure;
% view(3)
% hold on
% xlabel('$X$')
% ylabel('$Y$')
% zlabel('$Z$')
% plot3(0,0,0, '*k');
% plot3(x(1), y(1), z(1), '*k');
% % plot3(xE, yE, zE,'LineStyle','--','Color','r','LineWidth',0.3);   % Initial orbit
% % plot3(xM, yM, zM,'LineStyle','-.','Color','b','LineWidth',0.3);   % Final orbit
% hold on
% grid on; 
%     
% legend('off')
% plot3(x, y, z,'k','LineWidth',1);
% plot3(x(end), y(end), z(end),'*k');
% grid on;

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

figure 
hold on
plot(tau, I); 
hold off 
grid on;
xlabel('$t$')
ylabel('$A_i$')
title('Evolution of the integrals of motion')

figure 
hold on
plot(I(1,:), I(2,:)); 
hold off 
grid on;
xlabel('$A_x$')
ylabel('$A_z$')
title('Configuration space')

%% Auxiliary functions 
% Compute the Legendre coefficients 
function [c] = legendre_coefficients(mu, L, gamma, order)
    
    % Main computation 
    switch (L)
        case 1
            deg = 2:order; 
            c = (1 / gamma^3) .* ( mu + (1-mu) * (-1).^deg .* gamma.^(deg+1) ./ (1-gamma).^(deg+1) );

        case 2
            deg = 2:order; 
            c = ( (-1).^deg ./ gamma^3 ) .* ( mu + (1-mu) * gamma.^(deg+1) ./ (1+gamma).^(deg+1) );

        case 3
            deg = 2:order; 
            c = ( (-1).^deg ./ gamma^3 ) .* ( 1 - mu + ( mu * gamma.^(deg+1) ./ (1+gamma).^(deg+1) ) );

        otherwise
            error('No valid Lagrange point was selected');
    end

    % Add the zero for easy indexing
    c = [0 c];
end

% Compute the location of the libration points
function [L] = libration_points(mu)
    % Preallocation 
    colL = zeros(4,3); 
    
    % Compute the equilateral points (forming two symmetric equilateral triangles with the primaries)
    alpha = pi/3;                               % Angle between the libration points and primaries             
    primR = [-mu -mu; 0 0; 0 0; 1 1];           % Position of the first primary
    equiL = primR+[cos(alpha) cos(alpha); 
                   sin(alpha) -sin(alpha); 
                                  0 0; 0 0];    % Libration points positions
    
    % Compute the collinear points using a standard Newton method
    numL = 2;           % Number of collinear points to calculate recursively
    tol = 1e-15;        % Newton method tolerance
    iterMax = 1e5;      % Maximum allowed iterations for the Newton method
    
    % Main loop
    for i = 1:numL
        % Set up the Newton loop for L1/L2
        rh = mu^(1/3);                              % Hill radius
        lambda = rh*(1+((-1)^i)*(rh/3 +rh^2/9));    % Initial guess 
        GoOn = true;                                % Convergence flag
        iter = 1;                                   % Initial iteration
        
        % Main computation
        while ((GoOn) && (iter < iterMax))
            % Newton algorithm
            f = lambda^5 +((-1)^i)*(3-mu)*lambda^4 +(3-2*mu)*lambda^3 ...
                -mu*lambda^2 +2*((-1)^(i+1))*mu*lambda-mu;
            df = 5*lambda^4 +((-1)^i)*4*(3-mu)*lambda^3 +3*(3-2*mu)*lambda^2 ...
                 -2*mu*lambda +2*((-1)^(i+1))*mu;
            dn = -f/df;

            % Newton update
            lambda = lambda +dn;
            
            % Check for convergence
            if (abs(dn) < tol)
                GoOn = false;
            else
                iter = iter+1; 
            end
        end
        
        % Save the converged collinear point position in an array
        colL(:,i) = [(1-mu)+((-1)^(i))*lambda(end); 0; 0; lambda];
    end
    
    % Set up the Newton loop for L3
    rh = 1-(7/12)*mu;       % Initial guess
    lambda = rh;            % Initial guess
    GoOn = true;            % Convergence flag
    iter = 1;               % Initial iteration
    
    % Main computation
    while ((GoOn) && (iter < iterMax))
        f = lambda^5 +(2+mu)*lambda^4 +(1+2*mu)*lambda^3 -(1-mu)*lambda^2 -2*(1-mu)*lambda -(1-mu);
        df = 5*lambda^4 +4*(2+mu)*lambda^3 +3*(1+2*mu)*lambda^2 -2*(1-mu)*lambda -2*(1-mu);

        % Newton update
        dn = -f/df;
        lambda = lambda +dn;

        % Convergence check
        if (abs(dn) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end
    colL(:,3) = [-(mu+lambda); 0; 0; lambda];
    
    % Save results in the L structure ouput 
    L = [colL equiL];
end