%% Project: SBOPT %%
% Date: 18/01/23

%% Close-range rendezvous in YA model %% 
% This script provides the solving of Problem IIb) in the URJC rendezvous laboratory,
% SciTech2024

%% Set up 
close all
clear

%% Numerical solver definition 
basis = 'Bernstein';                    % Polynomial basis to be use
time_distribution = 'Bernstein';        % Distribution of time intervals
baseline_flag = true;                  % Baseline (true) vs reduced solution (false)

if (baseline_flag)
    N = 10;                            % Polynomial order in the state vector expansion
    m = 100;                           % Number of sampling points
else
    N = 10;                            % Polynomial order in the state vector expansion
    m = 200;                           % Number of sampling points
end
 
solver = Solver(basis, N, time_distribution, m);

%% Problem definition 
% Mission constraints
TOF = 6 * 3600;                      % Maximum allowed phase time [s]

% Target's initial conditions
ST = [0.5 * ones(4,1); deg2rad(1); deg2rad(0); deg2rad(-1)];

% Chaser's conditions
SC = [zeros(3,1); 1; zeros(3,1)];   

% Physical parameters
Fmax = 20e0;                            % Maximum available torque [m/s^2]
It = diag([140 36.9 36.9]);             % Inertia tensor of the target [kg m^2]
Ic = diag([100 95.9 80]);               % Inertia tensor of the chaser [kg m^2]

theta_e = deg2rad(10);                  % Error angle [rad]
omega_max = deg2rad(100);                % Maximum chaser's angular velocity infinity norm

% Final boundary conditions
SF = [ST; 0];       % Final conditions

%% Problem parameters
% Angular problem data
params(1) = TOF;                        % Maximum TOF [s]
params(2) = theta_e;                    % Attitude angular error [rad]
params(3) = omega_max;                  % Maximum angular velocity [rad/s]
params(4:12) = reshape(Ic, 3, []);      % Chaser inertia tensor
params(13:21) = reshape(It, 3, []);     % Target inertia tensor
params(22:28) = reshape(ST, 1, []);     % Initial target conditions 
params(29:31) = [0 0 1];                % Robotic arm boresight direction in the chaser body frame  
params(32:34) = [0 0 1];                % Graspling fixture boresight direction in the target body frame 
params(35) = Fmax;                      % Torque envelope [N m] 

%% Create the problem
L = 2;                                  % Degree of the dynamics (maximum derivative order of the ODE system)
StateDimension = 4;                     % Dimension of the configuration vector. Note the difference with the state vector
ControlDimension = 3;                   % Dimension of the control vector

%% Optimization (NMPC-RTI)
% Setup
options = odeset('AbsTol', 1e-22, 'RelTol', 2.25e-14);  % Integration tolerances
Ts = 600;                                                % Sampling time

% Noise 
Sigma_o = 0.01 * eye(3);                                % Relative angular velocity covariance

GoOn = true;                                            % Convergence boolean
iter = 1;                                               % Initial iteration
maxIter = 1e6;                                          % Maximum number of iterations
elapsed_time = 0;                                       % Elapsed time
tof = TOF;                                              % Maximum phase time

St = zeros(maxIter,7);
C = zeros(maxIter,7);
St(1,:) = ST.';
C(1,:) = SC.';
U = [];

% Target and chaser attitude conditions with respect to the LVLH frame
y0 = [ST; SC];  

while (GoOn && iter < maxIter)
    % Optimization (feedback phase)
    elapsed_time = (iter-1) * Ts;

    omega_0 = 0.5 * QuaternionAlgebra.right_isoclinic([C(iter,5:7).'; 0]) * C(iter,1:4).';   % Quaternion kinematics
    S0 = [C(iter,1:4).'; omega_0];
    OptProblem = Problems.RobotRotoBerthing(S0, SF, L, StateDimension, ControlDimension, params);
    [~, ~, u, t0, tf, tau, exitflag, ~, P] = solver.solve(OptProblem);
     
    % Control trajectory
    if (norm(u(:,1), 'inf') > Fmax)
        u(norm(u(:,1), 'inf') > Fmax,1) = Fmax;
    end
    U = [U u(:,1)];     

    % Preparing phase
    solver.InitialGuessFlag = true; 
    solver.P0 = P;
    solver.maxIter = 1;
    
    % Plant dynamics 
    [~, s] = ode113(@(t,s)euler_dynamics(Ic, It, t, s, U(:,end)), [0 Ts], y0, options);  

    % Update initial conditions and state vector
    y0 = s(end,:);
    St(iter+1,:) = s(end,1:7);                         % Target attitude conditions
    C(iter+1,:) = s(end,8:14);                         % Chaser attitude conditions

    % Navigation system
    omega = mvnrnd(St(iter+1,5:7), Sigma_o, 1);                     % Noisy state vector
    params(22:28) = reshape([St(iter+1,1:4) omega], 1, []);         % Initial target conditions

    % Convergence 
    if (elapsed_time >= TOF)
        GoOn = false;
    else
        % Shrinking horizon
        if (elapsed_time > 0.8 * TOF)
            tof = tof - Ts; 
            params(1) = tof;
        end

        % Update the number of iterations
        iter = iter + 1;
    end
end

U = [U zeros(3,1)];     % Final control vector
t = Ts * (0:iter);      % Elapsed time vector

C = C(1:iter+1,:).';
St = St(1:iter+1,:).';
          
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

e1 = [ones(1,length(t)); zeros(2,length(t))]; 
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
function [dq] = euler_dynamics(I, It, t, s, u)
    c = 1e-3; 

    dq(1:4) = 0.5 * QuaternionAlgebra.right_isoclinic( [s(5:7); 0] ) * s(1:4) + c * (1 - s(1:4).' * s(1:4)) * s(1:4);
    dq(5:7) = It \ cross(It * s(5:7), s(5:7));

    dq(8:11) = 0.5 * QuaternionAlgebra.right_isoclinic( [s(12:14); 0] ) * s(8:11) + c * (1 - s(8:11).' * s(8:11)) * s(8:11);
    dq(12:14) = I \ ( u + cross(I * s(12:14), s(12:14)) );

    dq = dq.';
end