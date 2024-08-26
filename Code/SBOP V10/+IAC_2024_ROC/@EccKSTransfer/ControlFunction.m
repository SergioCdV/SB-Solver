%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Constants 
    mu = params(1); 

    % Sundman transformation 
    r = dot(s(1:4,:), s(1:4,:), 1);                % Radius of the trajectory

    % Compute the energy 
    [~, alpha] = LegoKS.OscEnergy(mu, s, "Ecc");

    % Linear terms of the equations of motion
    a = s(9:12,:);                                 % Inertial acceleration field

    % Compute the control vector as a dynamics residual
    u = 2 * mu .* alpha .* (a + 0.25 .* s(1:4,:));
    
    aux = 4 ./ (r + 4 * dot(s(5:8,:), s(5:8,:), 1));
    for i = 1:size(u,2)
        L = LegoKS.KSmatrix( s(1:4,i) );
%         B = eye(4) + 4 / r(i) * (s(5:8,i) * s(5:8,i).');
        Binv = eye(4) - aux(i) * (s(5:8,i) * s(5:8,i).');
        u(:,i) = L * ( Binv * u(:,i) );
    end
end