%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    omega = 2 * [s(4,:).*s(5,:)+s(3,:).*s(6,:)-s(2,:).*s(7,:)-s(1,:).*s(8,:); ...
            -s(3,:).*s(5,:)+s(4,:).*s(6,:)+s(1,:).*s(7,:)-s(2,:).*s(8,:); ...
             s(2,:).*s(5,:)-s(1,:).*s(6,:)+s(4,:).*s(7,:)-s(3,:).*s(8,:)];

    dq = QuaternionAlgebra.right_isoclinic(s(1:4,end)) * QuaternionAlgebra.quaternion_inverse( params(15:18).' );
    c = [
         dot(u(1:3,:), u(1:3,:),1)-params(3)^2 ...  % Constraint on the torque magnitude second order cone
         reshape(+omega-params(4), 1, []) ...       % Constraint on the maximum angular velocity (epigraph form)
         reshape(-omega-params(4), 1, []) ...       % Constraint on the maximum angular velocity (epigraph form)
         -dq(4) - cos(params(14)/2) ...             % Terminal attitude constraint
         +dq(4) - cos(params(14)/2);                % Terminal attitude constraint
         ];                                                                   

    % Equality constraints
    ceq = [
          (dot(s(1:4,:), s(1:4,:), 1).'-1).' ...    % Quaternion norm constraint
          omega(:,end).'-params(19:21) ...          % Terminal angular velocity 
          ];     
end