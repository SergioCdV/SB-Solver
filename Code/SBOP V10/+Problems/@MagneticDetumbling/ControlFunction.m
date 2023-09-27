%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Parameters
    I = reshape(params(5:13), [3 3]);
    u = zeros(3,size(s,2));

    % Compute the control vector as a dynamics residual
    for i = 1:size(s,2)
        u(:,i) = (params(1) * I * s(4:6,i) + cross( s(1:3,i), I * s(1:3,i) )) / beta(end);
    end

%     u = beta(end) * u;
end