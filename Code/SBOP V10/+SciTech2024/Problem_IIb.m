%% Project: SBOPT %%
% Date: 18/01/23

%% Close-range rendezvous in YA model %% 
% This script provides the solving of Problem IIb) in the URJC rendezvous laboratory,
% SciTech2024

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Legendre';                    % Polynomial basis to be use
time_distribution = 'Legendre';        % Distribution of time intervals
N = 10;                                % Polynomial order in the state vector expansion
m = 100;                               % Number of sampling points
 
solver = Solver(basis, N, time_distribution, m);

%% Problem definition 
% Mission constraints
TOF = 0.25 * 3600;                      % Maximum allowed phase time [s]

% Target's initial conditions
ST = [1; 0; 0; deg2rad(0); deg2rad(0); deg2rad(0)];

% Chaser's conditions
SC = [zeros(3,1); zeros(3,1)];   

% Physical parameters
Fmax = 1e1;                            % Maximum available torque [m/s^2]
It = diag([140 36.9 36.9]);             % Inertia tensor of the target [kg m^2]
Ic = diag([100 95.9 80]);               % Inertia tensor of the chaser [kg m^2]

theta_e = deg2rad(10);                  % Error angle [rad]
omega_max = deg2rad(10);                % Maximum chaser's angular velocity infinity norm

% Final boundary conditions
SF = ST;                           % Final conditions

%% Dimensionalisation
Ip = It(1,1);
Ic = Ic / Ip;
Fmax = Fmax / Ip;
It = It / Ip;

%% Problem parameters
% MPC data 
Th = 0.1 * TOF;                         % Prediction horizon [s]
Ts = 10;                          % Control horizon [s]

% Problem's data
params(1) = 0;                          % Initial clock [s]
params(2) = Th;                         % Maximum TOF [s]
params(3) = theta_e;                    % Attitude angular error [rad]
params(4) = omega_max;                  % Maximum angular velocity [rad/s]
params(5) = Fmax;                       % Torque envelope [N m]
params(6:14) = reshape(Ic, 3, []);      % Chaser inertia tensor
params(15:23) = reshape(It, 3, []);     % Target inertia tensor
params(24:29) = reshape(ST, 1, []);     % Initial target conditions 
params(30:32) = [0 0 1];                % Robotic arm boresight direction in the chaser body frame  
params(33:35) = [0 0 1];                % Graspling fixture boresight direction in the target body frame 
 
%% Create the problem
L = 2;                                  % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 3;                     % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;                   % Dimension of the control vector

%% Optimization (NMPC-RTI)
% Setup
options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);  % Integration tolerances

% Noise 
Sigma_c = (deg2rad(0.1))^2 * eye(3);                    % Chaser's angular velocity covariance
Sigma_o = (deg2rad(0.1))^2 * eye(3);                    % Target's angular velocity covariance

GoOn = true;                                            % Convergence boolean
iter = 1;                                               % Initial iteration
maxIter = 1e6;                                          % Maximum number of iterations
elapsed_time = 0;                                       % Elapsed time

St = zeros(maxIter,6);
C = zeros(maxIter,6);
St(1,:) = ST.';
C(1,:) = SC.';
U = [];

% Target and chaser attitude conditions with respect to the LVLH frame
y0 = [ST; SC];  

while (GoOn && iter < maxIter)
    % Optimization (feedback phase)
    % omega_0 = 0.5 * QuaternionAlgebra.right_isoclinic([C(iter,5:7).'; 0]) * C(iter,1:4).';   % Quaternion kinematics
    S0 = C(iter,:).';

    if (norm(S0(1:3)) > 1)
        S0(1:3) = - S0(1:3) / dot(S0(1:3), S0(1:3));
    end

    tic
    OptProblem = Problems.RobotRotoBerthing(S0, SF, L, StateDimension, ControlDimension, params);
    [T, dV, u, t0, tf, tau, exitflag, ~, P] = solver.solve(OptProblem);

    if (iter == 1)
        %   params(36:43) = T(1:8,end).';
        open_time = toc();
        close_cost = dV;
    else
        % Save the computational delay
        iter_time(iter-1) = toc();
        
        % Integration during the preparation phase
        [~, s] = ode45(@(t,s)euler_dynamics(Ic, It, t, s, Pu, t0, tf), [Ts Ts + iter_time(iter-1)], y0, options);
        U = [U Pu * PolynomialBases.Legendre().basis(size(Pu,2)-1, 2* (linspace(Ts, Ts + iter_time(iter-1), 100) - t0) / (tf-t0) -1)];
        y0 = s(end,:);
        elapsed_time = elapsed_time + iter_time(iter-1);
    end
     
    % Control trajectory
    u(:, sqrt(dot(u,u,1)) > Fmax) = Fmax;
    Pu = PolynomialBases.Legendre().modal_projection(u);
    U = [U Pu * PolynomialBases.Legendre().basis(size(Pu,2)-1, 2* (linspace(0, Ts, 100) - t0)/(tf-t0) -1)];     

    % Preparing phase
    solver.InitialGuessFlag = true; 
    solver.P0 = P;
    solver.maxIter = 10;
    
    % Plant dynamics 
    [~, s] = ode45(@(t,s)euler_dynamics(Ic, It, t, s, Pu, t0, tf), [0 Ts], y0, options); 
    elapsed_time = elapsed_time + Ts;

    % Update initial conditions and state vector
    y0 = s(end,:);                 % New initial conditions
    St(iter+1,:) = s(end,1:6);     % Target attitude conditions
    C(iter+1,:) = s(end,7:12);     % Chaser attitude conditions
    
    if (norm(y0(1:3)) > 1)
        y0(1:3) = -y0(1:3) / dot(y0(1:3),y0(1:3));
    end

    if (norm(y0(7:9)) > 1)
        y0(7:9) = -y0(7:9) / dot(y0(7:9),y0(7:9));
    end

    % Navigation system
    omega = mvnrnd(St(iter+1,4:6), Sigma_o, 1);                     % Noisy state vector
%     params(22:28) = reshape([St(iter+1,13) omega], 1, []);         % Initial target conditions
%     C(iter+1,4:6) = mvnrnd(C(iter+1,4:6), Sigma_c, 1);              % Initial chaser conditions

    % Convergence 
    if (elapsed_time >= TOF)
        GoOn = false;
    else
        % Update the number of iterations
        iter = iter + 1;
    end
end

U = [U zeros(3,1)];     % Final control vector
t = Ts * (0:iter);      % Elapsed time vector

C = C(1:iter+1,:).';
St = St(1:iter+1,:).';

iter_time = mean(iter_time);
          
%% Plots
% Relative state representation
qr = C(1:4,:); 
omega_r = C(5:7,:); 
for i = 1:length(t)
    aux = QuaternionAlgebra.RotateVector(qr(:,i), St(5:7,i));
    omega_r(:,i) = C(5:7,i) - aux(1:3);

    qr(:,i) = QuaternionAlgebra.right_isoclinic(C(1:4,i)) * QuaternionAlgebra.quaternion_inverse(St(1:4,i));
end

figure
subplot(1,2,1)
hold on
xlabel('$t$ [min]')
ylabel('$\mathbf{q}$', 'Interpreter', 'latex')
plot(t / 60, qr(1:4,:));
legend('$q_1$', '$q_2$', '$q_3$', '$q_4$')
hold off
grid on;
xlim([0 t(end) / 60])

subplot(1,2,2)
hold on
xlabel('$t$ [min]')
ylabel('$\dot{\mathbf{\omega}}$ [rad/s]')
plot(t / 60, omega_r );
legend('$\dot{\omega}_x$', '$\dot{\omega}_y$', '$\dot{\omega}_z$')
hold off
grid on;
xlim([0 t(end) / 60])

% Propulsive acceleration plot
figure;
hold on
plot(t / 60, U(1:3,:), 'LineWidth', 0.3)
plot(t / 00, max(abs(U(1:3,:)), [], 1), 'k');
yline(Fmax, 'k--')
xlabel('$t$ [min]')
ylabel('$\mathbf{u}$ [N m]')
legend('$\tau_x$', '$\tau_y$', '$\tau$', '$\|\mathbf{\tau}\|_\infty$', '$\tau_{max}$');
grid on;
xlim([0 t(end)] / 60)

%% Sphere figure 
figure 
view(3)
hold on

no = [eye(3); zeros(1,3)];
nref = eye(3);
for j = 1:size(no,2)
    nref(:,j) = QuaternionAlgebra.RotateVector(qr(:,1), no(:,j));
end

quiver3(zeros(1,size(no,2)), zeros(1,size(no,2)), zeros(1,size(no,2)), no(1,:), no(2,:), no(3,:), 'b', 'LineWidth', 1.0);
quiver3(zeros(1,size(nref,2)), zeros(1,size(nref,2)), zeros(1,size(nref,2)), nref(1,:), nref(2,:), nref(3,:), 'r', 'LineWidth', 1.0);

for j = 1:size(no,2)
    nref(:,j) = QuaternionAlgebra.RotateVector(qr(:,end), no(:,j));
end

quiver3(zeros(1,size(nref,2)), zeros(1,size(nref,2)), zeros(1,size(nref,2)), nref(1,:), nref(2,:), nref(3,:), 'g', 'LineWidth', 1.0);

e1 = [zeros(2,length(t)); ones(1,length(t))]; 
for i = 1:length(t)
    e1(:,i) = QuaternionAlgebra.RotateVector(qr(:,i), e1(:,1));
end
legend('$\mathcal{B}_T$', '$\mathcal{B}_C(0)$', '$\mathcal{B}_C(t_f)$', 'AutoUpdate', 'off', 'Location', 'best')

scatter3(e1(1,:), e1(2,:), e1(3,:), 'filled')

m = 10;
[aa, bb, cc] = sphere(m);
h = surf(aa, bb, cc);
set(h, 'FaceColor', [0 0 1], 'FaceAlpha', 0.02, 'FaceLighting', 'gouraud', 'EdgeColor', 'none')
hold off
xlabel('$x_T$')
ylabel('$y_T$')
zlabel('$z_T$')
grid on;
xticklabels(strrep(xticklabels, '-', '$-$'));
yticklabels(strrep(yticklabels, '-', '$-$'));
zticklabels(strrep(zticklabels, '-', '$-$'));

%% Auxiliary functions 
% Attitude dynamics
function [dq] = euler_dynamics(I, It, t, s, P, t0, tf)
%     c = 0 * 1e-3; 

    tau = 2 * (t-t0)/(tf-t0) - 1;
    B = PolynomialBases.Legendre().basis(size(P,2)-1, tau);
    u = P * B;

%     dq(1:3) = 0.25 * QuaternionAlgebra.right_isoclinic( [s(5:7); 0] ) * s(1:4) + c * (1 - s(1:4).' * s(1:4)) * s(1:4);
%     dq(7:9) = 0.25 * QuaternionAlgebra.right_isoclinic( [s(12:14); 0] ) * s(8:11) + c * (1 - s(8:11).' * s(8:11)) * s(8:11);

    q1 = [s(1:3); -1];
    q2 = [s(7:9); -1];
    
    B1 = QuaternionAlgebra.Quat2Matrix(q1);
    B2 = QuaternionAlgebra.Quat2Matrix(q2); 

    dq(1:3) = 0.25 * B1 * s(4:6);
    dq(7:9) = 0.25 * B2 * s(10:12);

    dq(4:6) = It \ cross(It * s(4:6), s(4:6));
    dq(10:12) = I * ( u + cross(I * s(10:12), s(10:12)) );

    dq = dq.';
end