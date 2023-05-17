%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Constants 
    mu = params(1); 
    
    % Compute the radius vector
    r = s(1,:);

    % Linear terms of the equations of motion
    f = -[mu.*s(1,:)./r.^3; zeros(1,size(s,2))];                            % Acceleration vector
    a = [s(5,:)-s(1,:).*s(4,:).^2; s(1,:).*s(6,:)+2*s(3,:).*s(4,:)];        % Inertial acceleration field

    % Compute the control vector as a dynamics residual
    u = a-f;
end