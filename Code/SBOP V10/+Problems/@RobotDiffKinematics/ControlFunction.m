%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Preallocation 
    u = zeros(6, size(tau,2));

    % Evaluate the reference control law
    s_ref = reshape(params(43:end), [], size(tau,2));

    % Compute the Jacobian
    for i = 1:length(tau)
        [~, J] = Problems.RobotDiffKinematics.Kinematics(obj.StateDim, ...
                                                        @(j,s)Problems.RobotDiffKinematics.ur3_dkinematics(obj, j, s), ...
                                                        s(:,i));

        det(J)
        % Jacobian pseudoinverse
        invJ = pinv(J);
    end

    % Compute the control vector
    u(:,i) = s(7:12,i);
end