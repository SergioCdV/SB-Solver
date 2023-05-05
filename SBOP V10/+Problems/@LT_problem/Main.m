%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/08/22

%% Set up
set_graphics(); 
close all

%% Problem definition 
% Numerical solver definition 
time_distribution = 'Bernstein';       % Distribution of time intervals. 
basis = 'Bernstein';                   % Polynomial basis to be use
n = 7;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
L = 2;                                 % Degree of the dynamics 

n = 15; 

OptProblem = Problem().DefineSolver(n, basis, m, time_distribution).AddDynamics(1, 1, L); 

% Boundary conditions
S0 = [100; 20];
SF = [0; 0];

T = 1;              % Maximum acceleration 
g = 4;              % Acceleration

% Add boundary conditions
OptProblem = OptProblem.AddBoundaryConditions(S0, SF).AddParameters([g; T]);

% Add functions 
OptProblem = OptProblem.AddFunctions(@(initial, final, beta, t0, tf)BoundaryConditionsML(initial, final, beta, t0, tf), @(params, beta, t0, tf, tau, s)ControlFunctionML(params, beta, t0, tf, tau, s), ...
                                     @(params, beta, t0, tf, s, u)CostFunctionML(params, beta, t0, tf, s, u), @(beta, P)LinConstraintsML(beta, P), ...
                                     @(params, beta, t0, tf, tau, s, u)NlinConstraintsML(params, beta, t0, tf, tau, s, u), ...
                                     @BoundsFunctionML, ...
                                     @(params, initial, final)InitialGuessML(params, initial, final));

%% Optimization
% Simple solution    
tic
[C, dV, u, t0, tf, tau, exitflag, output] = sb_solver(OptProblem);
toc 

dV = dV * Vc;

% Average results 
iter = 0; 
time = zeros(1,iter);
setup.resultsFlag = false; 
for i = 1:iter
    tic 
    [C, dV, u, t0, tf, tau, exitflag, output] = sb_solver(OptProblem);
    time(i) = toc;
end

time = mean(time);

%% Plots
% Main plots 
[S] = cylindrical2cartesian(C(1:3,:),true);
x = S(1,:);
y = S(2,:); 
z = S(3,:);
    
% Earth's orbit
thetaE = linspace(0, 2*pi, size(C,2));

s = coe2state(mu, initial_coe);
initial = cylindrical2cartesian(s, false).';

s = coe2state(mu, final_coe);
final = cylindrical2cartesian(s, false).';
    
s = zeros(6,length(thetaE));
for i = 1:length(thetaE)
    s(:,i) = coe2state(mu, [initial_coe(1:end-1) initial(2)+thetaE(i)]);
end
xE = s(1,:);
yE = s(2,:);
zE = s(3,:);
    
% Mars's orbit
for i = 1:length(thetaE)
    s(:,i) = coe2state(mu, [final_coe(1:end-1) final(2)+thetaE(i)]);
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

%% Auxiliary functions
% Cylindrical to Cartesian
function [S] = cylindrical2cartesian(s, direction)
    % Sanity check on the s dimensions 
    if (size(s,1) == 3)
        lastwarn('State vector has only 3 dimensions')
        s = [s; zeros(3,size(s,2))];
    end

    % Switch directions 
    if (direction)
        % Cylindrical position coordinates
        rho = s(1,:);
        theta = s(2,:);
        z = s(3,:);

        % Cylindrical velocity coordinates
        drho = s(4,:);
        dtheta = s(5,:);
        dz = s(6,:);

        % Cartesian position coordinates
        S(1,:) = rho.*cos(theta);
        S(2,:) = rho.*sin(theta);
        S(3,:) = z;

        % Cartesian velocity coordinates 
        S(4,:) = drho-rho.*sin(theta).*dtheta;
        S(5,:) = drho+rho.*cos(theta).*dtheta;
        S(6,:) = dz;
    
    else
        % Cartesian position coordinates
        x = s(1,:);
        y = s(2,:);
        z = s(3,:);

        % Cartesian velocity coordinates
        dx = s(4,:);
        dy = s(5,:);
        dz = s(6,:);

        % Cartesian position coordinates
        S(1,:) = sqrt(x.^2+y.^2);
        S(2,:) = atan2(y,x);
        S(3,:) = z;

        % Cartesian velocity coordinates 
        S(4,:) = (x.*dx+y.*dy)./S(1,:);
        S(5,:) = (x.*dy-y.*dx)./S(1,:).^2;
        S(6,:) = dz;
    end  
end

% State to Cartesian 
function [elements] = state2coe(mu, s, frame)
    % Main computation 
    switch (frame)
        case 'Inertial'
            elements = rv2coe(mu,s);
        case 'Perifocal' 
            elements  = prv2coe(mu,s);
        otherwise
            error('No valid reference frame is used');
    end
end

function [s] = coe2state(mu, elements)
    % Constants 
    e = elements(2);                                            % Eccentricity of the orbit
    
    % Singularity warnings 
    tol = 1e-10;                                                % Circular orbit tolerance    
    if (abs(norm(e)) < tol)
        elements(5) = 0;
    end
    
    if (elements(4) == 0)
        elements(3) = 0;
    end
    
    % Compute the semilatus rectum
    if (elements(1) == Inf) 
        p = elements(end);                                      % Semilatus rectum of the orbit
    else
        p = elements(1)*(1-elements(2)^2);                      % Semilatus rectum of the orbit
    end
    
    % Compute the angular momentum norm
    h = sqrt(mu*p);                                             % Angular momentum of the orbit
    
    % Compute the mean anomaly
    theta = kepler(elements);                                   % True anomaly in the orbit
    
    % Compute the perifocal state vector
    r = (p/(1+e*cos(theta)))*[cos(theta); sin(theta); 0];       % Position vector in the perifocal frame
    v = mu/h*[-sin(theta); e+cos(theta); 0];                    % Velocity vector in the perifocal frame
    
    % Rotation matrix from the inertial to the perifocal frame
    Q = euler_matrix(elements);
       
    % Output
    r = Q.'*r;      % Position vector in the inertial frame
    v = Q.'*v;      % Velocity vector in the inertial frame
    s = [r; v];     % State vector
end
 
% Function to compute the orbital elements from the inertial state vector 
function [elements] = rv2coe(mu, s)
    % State variables 
    r = s(1:3);                                     % Position vector
    v = s(4:6);                                     % Velocity vector 
    
    % Compute the eccentricy and angular momentum vectors 
    h = cross(r,v);                                 % Angular momentum vector
    e = cross(v,h)/mu-r/norm(r);                    % Eccentricity vector
    K = [0; 0; 1];                                  % Inertial Z axis unit vector
    n = cross(K, h);                                % Node vector
    
    % Compute orbital energy 
    H = norm(v)^2/2-mu/norm(r);
    
    % Determine type of orbit 
    if (norm(e) ~= 1)
        a = -mu/(2*H);                              % Semimajor axis of the orbit
        p = a*(1-norm(e)^2);                        % Semilatus rectum of the orbit
    else
        p = norm(h)^2/mu;                           % Semilatus rectum of the orbit
        a = Inf;                                    % Semimajor axis of the orbit
    end
    
    % Compute the unit perifocal triad 
    i = e/norm(e); 
    k = h/norm(h); 
    j = cross(k,i);
    
    % Compute the rotation matrix
    Q = [i.'; j.'; k.'];                            % Perifocal rotation matrix
    r0 = Q*r;                                       % Position in the perifocal frame 
    
    % Compute the rest of the orbital elements
    RAAN = atan2(Q(3,1),-Q(3,2));                                        % RAAN
    omega = atan2(Q(1,3),Q(2,3));                                        % Argument of perigee
    i = acos(Q(3,3));                                                    % Inclination of the orbit
    
    % Mean anomaly
    theta = atan2(r0(2), r0(1));                                         % True anomaly of the orbit
    sinE = sqrt(1-norm(e)^2)*sin(theta)/(1+norm(e)*cos(theta));          % Sine of the eccentric anomaly
    cosE = (norm(e)+cos(theta))/(1+norm(e)*cos(theta));                  % Cosine of the eccentric anomaly
    E = atan2(sinE, cosE);                                               % Eccentric anomaly
    M = E-norm(e)*sin(E);                                                % Mean anomaly
        
    % Save the classical orbital elements 
    elements = [a norm(e) RAAN i omega M p];
    elements = rv_singularity(e, n, r, Q, elements);
end

% Function to compute the orbital elements from the perifocal state vector
function [elements] = prv2coe(mu, s)
    % State variables 
    r = s(1:3).';                                   % Position vector
    v = s(4:6).';                                   % Velocity vector 
    
    % Compute the eccentricy and angular momentum vectors 
    h = cross(r,v);                                 % Angular momentum vector
    h = norm(h);                                    % Angular momentum norm
    
    % Compute the true anomaly 
    theta = atan2(r(2), r(1)); 
    
    % Compute the eccentricity 
    e = v(2)*(h/mu)-cos(theta); 
    
    % Compute orbital energy 
    H = norm(v)^2/2-mu/norm(r);
    
    % Determine type of orbit 
    if (e ~= 1)
        a = -mu/(2*H);                              % Semimajor axis of the orbit
        p = a*(1-e^2);                              % Semilatus rectum of the orbit
    else
        p = h^2/mu;                                 % Semilatus rectum of the orbit
        a = Inf;                                    % Semimajor axis of the orbit
    end
        
    % Compute the rest of the orbital elements
    RAAN = NaN;                                         % RAAN
    omega = NaN;                                        % Argument of perigee
    i = NaN;                                            % Inclination of the orbit
    
    % Mean anomaly
    sinE = sqrt(1-e^2)*sin(theta)/(1+e*cos(theta));     % Sine of the eccentric anomaly
    cosE = (e+cos(theta))/(1+e*cos(theta));             % Cosine of the eccentric anomaly
    E = atan2(sinE, cosE);                              % Eccentric anomaly
    M = E-norm(e)*sin(E);                               % Mean anomaly
    
    % Save the classical orbital elements 
    elements = [a e RAAN i omega M p];
    
    % Singularity warnings 
    tol = 1e-10;               % Circular orbit tolerance    
    if (abs(e) < tol)
        warning('Orbit is circular to numerical precision');   
    end
end

% Function to handle orbital elements singularities when converting from the inertial state vector
function [elements] = rv_singularity(ev, n, r, Q, elements)
    % State variables 
    e = elements(2);           % Orbit eccentricity
    i = elements(4);           % Orbit inclination 
    M = elements(6);           % Mean anomaly
    tol = 1e-10;               % Circular orbit tolerance
    
    % Singularity warnings 
    if (any(isnan(Q)))
        warning('Euler angles are numerically ill-conditioned');
    end
    
    if (abs(e) < tol)
        warning('Orbit is circular to numerical precision');
    end
    
    if (abs(i) < tol)
        warning('Orbit is equatorial to numerical precision');
    end
    
    %Singularity handling 
    if (abs(e) < tol)
        if (abs(i) < tol)
            M = acos(r(1)/norm(r));                 % Circular equatorial orbit
            if (r(2) < 0)
                M = 2*pi-M;
            end
        else
            M = acos(dot(n,r)/(norm(n)*norm(r)));   % Circular inclined orbit
            if (r(3) < 0)
                M = 2*pi-M;
            end
        end
    else
        if (abs(i) < tol)
            M = acos(ev(1)/norm(ev));               % Equatorial orbit
            if (ev(2) < 0)
                M = 2*pi-M; 
            end
        end
    end
    
    % Reconstruction of the orbital elements 
    elements(6) = M;
end

% Kepler equation 
function [theta] = kepler(elements)
    % Constants 
    e = elements(2);                                            % Eccentricity of the orbit
    M0 = elements(6);                                           % Mean anomaly of the orbit 
    
    % Set up the loop 
    k = 5;                                                      % Conway constant
    tol = 1e-15;                                                % Convergence tolerance
    iterMax = 10^6;                                             % Maximum number of iterations
    GoOn = true;                                                % Convergence flag
    iter = 1;                                                   % Initial iteration
    u = M0+e;                                                   % Conway method variable
    E(iter) = (M0*(1-sin(u))+u*sin(M0))/(1+sin(M0)-sin(u));     % Initial guess for the eccentric anomaly
    
    % Main computation 
    while ((GoOn) && (iter < iterMax))
        % Laguerre-Conway iterations
        f = E(iter)-e*sin(E(iter))-M0; 
        df = 1-e*cos(E(iter));
        ddf = e*sin(E(iter));
        dn = -k*(f)/(df+sqrt((k-1)^2*df^2-k*(k-1)*f*ddf));
        E(iter+1) = E(iter)+dn;
        
        % Convergence checking 
        if (abs(dn) < tol)
            GoOn = false;
        else
            iter = iter+1;
        end
    end  
    
    % True anomaly
    theta = atan2(sqrt(1-e^2)*sin(E(end))/(1-e*cos(E(end))), (cos(E(end))-e)/(1-e*cos(E(end))));
end

% Euler matrix 
function [Q] = euler_matrix(elements)
    % Elements of interest 
    RAAN = elements(3); 
    i = elements(4); 
    omega = elements(5); 
    
    % Compute the rotation matrix (Euler sequence ZXZ)
    Q1 = [cos(RAAN) sin(RAAN) 0; -sin(RAAN) cos(RAAN) 0; 0 0 1];
    Q2 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    Q3 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
    Q = Q3*Q2*Q1;
end