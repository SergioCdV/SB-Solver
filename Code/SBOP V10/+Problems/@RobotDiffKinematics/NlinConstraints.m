%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    c = [-reshape(u, 1, [])-params(3) ... % Constraint on the torque magnitude (infinity norm in epigraph form)
          reshape(u, 1, [])-params(3) ... % Constraint on the torque magnitude (infinity norm in epigraph form)
        ];                                                                                      

    % Equality constraints
    base = params(4:6).'; 
    theta = params(7:12).';
    alpha = params(13:18).';
    offset = params(19:24).';
    a = params(25:30).';
    d = params(31:36).';
    type = params(37:42).';        % All joints are revolute

    [T, J] = Problems.RobotDiffKinematics.Kinematics(obj.StateDim, type, ...
                                                     @(i,s)Problems.RobotDiffKinematics.ur3_dkinematics(obj, base, theta, alpha, offset, a, d, type, i, s), ...
                                                     s(:,end));
   
    s_ref = reshape(params(43:end), [], size(tau,2));   % Final end-effect operational state
    q = T * [0; 0; 1; 1];                               % End-effector position
    v = J * s(obj.StateDim+1:2*obj.StateDim,end);       % End-effector linear and angular velocity

    bcs = [[0; 0; 1; 1]-pinv(T) * [s_ref(1:3,end); 1];...
           [v(1:3); v(4:6)]-s_ref([8:end],end) ...
          ];       % Boundary conditions constraint

    ceq = bcs;    
end