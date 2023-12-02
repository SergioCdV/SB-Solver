%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, tau, s)
    % Compute the angular velocity
    I = reshape(params(6:14), 3, 3);      % Inertia tensor of the chaser

%     % Angular velocity
%     omega = 2 * [s(4,:).*s(5,:)+s(3,:).*s(6,:)-s(2,:).*s(7,:)-s(1,:).*s(8,:); ...
%                 -s(3,:).*s(5,:)+s(4,:).*s(6,:)+s(1,:).*s(7,:)-s(2,:).*s(8,:); ...
%                  s(2,:).*s(5,:)-s(1,:).*s(6,:)+s(4,:).*s(7,:)-s(3,:).*s(8,:)];
% 
%     % Angular acceleration 
%     alpha = 2 * [s(4,:).*s(9,:)+s(3,:).*s(10,:)-s(2,:).*s(11,:)-s(1,:).*s(12,:); ...
%                 -s(3,:).*s(9,:)+s(4,:).*s(10,:)+s(1,:).*s(11,:)-s(2,:).*s(12,:); ...
%                  s(2,:).*s(9,:)-s(1,:).*s(10,:)+s(4,:).*s(11,:)-s(3,:).*s(12,:)];

    omega = zeros(3,length(tau));
    alpha = omega;

    % Euler equations
    u = zeros(3,size(tau,2));
    for i = 1:length(tau)
%         norm = dot(s(1:3,i),s(1:3,i));
%         if norm > 1.05 
%             s(1:3,i) = -s(1:3,i) / norm;
%         end
        q = [s(1:3,i); -1];
        B = QuaternionAlgebra.Quat2Matrix(q);
        omega(:,i) = 4 * B.' * s(4:6,i) / dot(q,q)^2;
        dB = dBmatrix(s(1:3,i), s(4:6,i));
        alpha(:,i) = 4 * B.' * (s(7:9,i) - 0.25 * dB * omega(:,i)) / dot(q,q)^2;
        u(1:3,i) = I * alpha(:,i) + cross( omega(:,i), I*omega(:,i) );           
    end
end

%% Auxiliary function 
function [Bdot] = dBmatrix(q, dq)
    Bdot = zeros(3,3);
    s = -2 * dot(q, dq);
    Bdot(1,1) = s + 4 * (q(1) * dq(1));
    Bdot(1,2) = 2 * (-dq(3) + q(1) * dq(2) + dq(1) * q(2));
    Bdot(1,3) = 2 * ( dq(2) + q(1) * dq(3) + dq(1) * q(3));
    Bdot(2,1) = 2 * ( dq(3) + q(1) * dq(2) + dq(1) * q(2));
    Bdot(2,2) = s + 4 * (q(2) * dq(2));
    Bdot(2,3) = 2 * (-dq(1) + q(2) * dq(3) + dq(2) * q(3));
    Bdot(3,1) = 2 * (-dq(2) + q(1) * dq(3) + dq(1) * q(3));
    Bdot(3,2) = 2 * ( dq(1) + q(2) * dq(3) + dq(2) * q(3));
    Bdot(3,3) = s + 4 * (q(3) * dq(3));
end