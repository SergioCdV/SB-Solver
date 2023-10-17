%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Preallocation 
    u = zeros(6, size(tau,2));

    % Evaluate the reference control law
    s_ref = reshape(params(43:end), [], size(tau,2));

    base = params(4:6).'; 
    theta = params(7:12).';
    alpha = params(13:18).';
    offset = params(19:24).';
    a = params(25:30).';
    d = params(31:36).';
    type = params(37:42).';        % All joints are revolute

    % Compute the Jacobian
    for i = 1:length(tau)
        [~, J] = Problems.RobotDiffKinematics.Kinematics(obj.StateDim, type, ...
                                                        @(i,s)Problems.RobotDiffKinematics.ur3_dkinematics(obj, base, theta, alpha, offset, a, d, type, i, s), ...
                                                        s(:,i));

        % Jacobian pseudoinverse
        invJ = pinv(J);
    
        % Compute the control vector
        u(:,i) = invJ * s_ref(8:13,i);
    end
end