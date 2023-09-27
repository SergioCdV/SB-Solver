%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, S)
    % Compute the auxiliary variables
    w = 1+S(2,:).*cos(S(6,:))+S(3,:).*sin(S(6,:));
    s = 1+S(4,:).^2+S(5,:).^2;

    % Linear terms of the equations of motion
    f = tf*[zeros(5,size(S,2)); sqrt(mu.*S(1,:)).*(w./S(1,:)).^2];                  % Acceleration vector
    a = S(7:12,:);                                                                  % Inertial acceleration field
    dv = S(7:12,:);                                                                 % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = zeros(3,size(S,2));

    % Tangential component
    alpha = 2*S(1,:)./w.*sqrt(S(1,:)/mu); 
    u(2,:) = (a(1,:)-f(1,:))./alpha;

    % Normal component
    beta = sqrt(S(1,:)/mu).*S./(2*w);
    u(3,:) = sqrt(a(4,:).^2+a(5,:).^2)./beta;
    u(3,:) = u(3,:).*sign(a(4,:))./sign(cos(S(6,:)));

    % Radial component
    delta = sqrt(S(1,:)/mu);
    for i = 1:size(S,2)
        B = OrbitalDynamics.MEE_matrix(mu, S(:,i));
        u(1,i) = (a(2,i)-B(2,2:3)*u(2:3,i))^2+(a(3,i)-B(3,2:3)*u(2:3,i))^2;
    end
    u(1,:) = sqrt(u(1,:))./delta;
end