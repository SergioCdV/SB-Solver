%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Constants 
    mu = params(1); 

    % Sundman transformation 
    r = dot(s(1:4,:), s(1:4,:), 1);     % Radius of the trajectory

    % Compute the energy 
    E = (4 * dot(s(5:8,:), s(5:8,:), 1) - 2 * mu) ./ r;
    
    % Linear terms of the equations of motion
    a = s(9:12,:);                                 % Inertial acceleration field

    % Compute the control vector as a dynamics residual
    u = (2 * a - E .* s(1:4,:));

    for i = 1:size(u,2)
        L = Problems.KSTransfer.KS_matrix(s(1:4,i));
        u(:,i) = L * u(:,i);
    end
end