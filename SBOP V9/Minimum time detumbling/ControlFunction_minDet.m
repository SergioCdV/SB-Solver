%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(params, beta, t0, tf, tau, s)
    % Compute the control vector as a dynamics residual
    I = reshape(params(3:end), [3 3]);
    u = zeros(3,size(s,2));
    for i = 1:size(s,2)
        u(:,i) = I * s(4:6,i) + cross( s(1:3,i),I * s(1:3,i) );
    end
end