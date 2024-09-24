%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Constants 
    mu = params(1);                                 % Gravitational constant of the system

    % Compute the energy 
    E = LegoKS.OscEnergy(mu, s);

    % Linear terms of the equations of motion
    a = s(9:12,:);                                 % Inertial acceleration field

    % Compute the control vector as a dynamics residual
    u = 2 * a - E .* s(1:4,:);

    for i = 1:size(u,2)
        L = LegoKS.KSmatrix( s(1:4,i) );
        u(:,i) = L * u(:,i);
    end
end