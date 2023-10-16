%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Preallocation 
    u = zeros(6, size(tau,2));

    % Evaluate the reference control law
    s_ref = reshape(params(4:end), [], size(tau,2));

    % Compute the Jacobian
    for i = 1:length(tau)
        J = Kinematics(obj.StateDim, params(4:4+obj.StateDim-1), @(i,s)UR13_dkinematics(i,s), s(:,i));
        invJ = pinv(J);
%         A = eye(size(J,2)) - invJ * J;
    
        % Compute the control vector
        u(:,i) = s(7:12,i) - invJ * s_ref(7:12,i);
    end
end