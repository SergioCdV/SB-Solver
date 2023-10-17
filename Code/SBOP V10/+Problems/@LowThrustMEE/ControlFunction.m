%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, S)
    % Constants
    mu = params(1);     % Gravitational parameter
    maxIter = 100;      % Maximum number of iterations 
    iter = 1;           % Initial iteration
    GoOn = true;        % Convergence variable

    % Preallocation 
    u = zeros(3,size(S,2));
    u_ref = 1e3 * ones(size(u));

    % Compute the auxiliary variables
    while (GoOn)
        w = 1+S(2,:).*cos(t(1,:))+S(3,:).*sin(t(1,:));
        s = 1+S(4,:).^2+S(5,:).^2;
        gamma = sqrt(S(1,:) * mu) .* (w./S(1,:)).^2 + sqrt(S(1,:)/mu)./w.*(S(4,:).*sin(t(1,:)) + S(5,:).*cos(t(1,:))) .* u(3,:);
    
        % Linear terms of the equations of motion
        f = zeros(5,size(S,2));                     % Acceleration vector
        a = gamma .* S(6:10,:);                     % Inertial acceleration field
    
        % Tangential component
        alpha = 2*S(1,:)./w.*sqrt(S(1,:)/mu); 
        u(2,:) = (a(1,:)-f(1,:))./alpha;
    
        % Normal component
        beta = sqrt(S(1,:)/mu).*s./(2*w);
        u(3,:) = sqrt(a(4,:).^2+a(5,:).^2)./beta;
        u(3,:) = u(3,:).*sign( a(4,:)/cos(t(1,:)) );
    
        % Radial component
        delta = sqrt(S(1,:)/mu);
        for i = 1:size(S,2)
            B = OrbitalDynamics.MEE_matrix(mu, t(1,i), S(1:5,i));
            u(1,i) = (a(2,i)-B(2,2:3)*u(2:3,i))^2+(a(3,i)-B(3,2:3)*u(2:3,i))^2;
        end
        u(1,:) = sqrt(u(1,:))./delta;

        GoOn = max(abs(u_ref(3,:)-u(3,:))) > params(5) && iter < maxIter;

        if (GoOn)
            u_ref = u;
            iter = iter + 1;
        end
    end
end