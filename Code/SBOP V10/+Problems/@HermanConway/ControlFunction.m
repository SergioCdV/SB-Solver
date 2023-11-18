%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Constants 
    mu = params(1); 
    
    % Compute the radius vector
    r = s(1,:);
    theta = s(2,:); 
    vr = s(3,:); 
    vt = s(4,:) .* s(1,:);

    % Compute the control vector as a dynamics residual
    f(1,:) = s(5,:) - vt.^2 ./ r + mu ./ r.^2;
    f(2,:) = r .* s(6,:) + vr .* s(4,:) + vr .* vt ./ r;
    u = f / params(3);
end