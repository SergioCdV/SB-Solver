%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Constants 
    mu = params(1);                                 % Gravitational parameter of the problem

    % Inequality constraints
    c = dot(u,u,1)-params(2)^2;                     % Thrust modulation

    % Monotony of the time function 
    w = 1+s(2,:).*cos(tau)+s(3,:).*sin(tau);
    res = sqrt(mu*s(1,:)).*(w./s(1,:)).^2;
    for i = 1:size(s,2)
        B = obj.control_input(mu, s(:,i)); 
        res(i) = 1/(res(i)+B(6,3)*u(3,i));
    end

    % Inequalities
    c = [c -res];

    % Equality constraints
    ceq = [cos(tf)-cos(params(3)); sin(tf)-sin(params(3))];     % Satisfaction of the boundary conditions

    % Inequalities
    c = [u(1,:).^2+u(2,:).^2+u(3,:).^2-(tf*repmat(T,1,size(u,2))).^2 -diff(C(6,:))];

    % Kinematic constraint 
    w = 1+C(2,:).*cos(C(6,:))+C(3,:).*sin(C(6,:));
    res = tf*sqrt(mu*C(1,:)).*(w./C(1,:)).^2;
    for i = 1:size(C,2)
        B = control_input(mu, C(:,i)); 
        res(i) = res(i)+B(6,3)*u(3,i);
    end


    % Equalities
    ceq = [cos(C(6,end))-cos(final(6)) sin(C(6,end))-sin(final(6))];
end