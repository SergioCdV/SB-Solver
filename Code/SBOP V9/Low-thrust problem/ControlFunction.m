%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(params, beta, t0, tf, t, s)
    % Constants 
    mu = params(1); 
    
    % Compute the radius vector
    r = sqrt(s(1,:).^2+s(3,:).^2);

    % Linear terms of the equations of motion
    c = (tf-t0);                                                                     % Normalizing factor
    f = -[mu.*s(1,:)./r.^3; zeros(1,size(s,2)); mu.*s(3,:)./r.^3];                  % Acceleration vector
    a = [s(7,:)-s(1,:).*s(5,:).^2; s(1,:).*s(8,:)+2*s(4,:).*s(5,:); s(9,:)];        % Inertial acceleration field

    % Compute the control vector as a dynamics residual
    u = a-c.^2.*f;
end