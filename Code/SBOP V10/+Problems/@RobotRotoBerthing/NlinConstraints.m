%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints    
    omega = 2 * [+s(4,:).*s(5,:)+s(3,:).*s(6,:)-s(2,:).*s(7,:)-s(1,:).*s(8,:); ...
                 -s(3,:).*s(5,:)+s(4,:).*s(6,:)+s(1,:).*s(7,:)-s(2,:).*s(8,:); ...
                 +s(2,:).*s(5,:)-s(1,:).*s(6,:)+s(4,:).*s(7,:)-s(3,:).*s(8,:)];

    % Relative state 
%     It = reshape(params(13:21), 3, 3);                  % Inertia tensor of the target
%     ST = reshape(params(22:28), [], 1);                 % Initial target state
%     Omega = (1 - It(1,1) / It(3,3)) * ST(7,1);
% 
%     omega_t = [repmat(ST(7,1), 1, length(tau)); ...
%                +ST(6) * cos(Omega * tau(1,:)) + ST(7) * sin(Omega * tau(1,:)); ...
%                -ST(6) * sin(Omega * tau(1,:)) + ST(7) * cos(Omega * tau(1,:))];
%     
%     % Target trajectory
%     q_t = s(1:4,:);
%     for i = 1:length(tau)
%         if (i == 1)
%             q_t(:,i) = ST(1:4,1);
%         else
%             q_t(:,i) = QuaternionAlgebra.right_isoclinic( QuaternionAlgebra.exp_map(0.5 * omega_t(:,i-1) * (tau(1,i)-tau(1,i-1)) ) ) * q_t(:,i-1);
%         end
%     end

    % Inequalities
%     r = params(29:31).';
%     b = params(32:34).'; 

%     r = QuaternionAlgebra.RotateVector( QuaternionAlgebra.quaternion_inverse(s(1:4,end)), r);
%     b = QuaternionAlgebra.RotateVector( QuaternionAlgebra.quaternion_inverse(q_t(:,end)), b);

    c = [
%              reshape(+u-params(35), 1, []) ...          % Constraint on the maximum control torque (epigraph form)
%              reshape(-u-params(35), 1, []) ...          % Constraint on the maximum control torque (epigraph form)
%              reshape(+omega-params(3), 1, []) ...       % Constraint on the maximum angular velocity (epigraph form)
%              reshape(-omega-params(3), 1, []) ...       % Constraint on the maximum angular velocity (epigraph form)
%                cos(params(2)) - dot(-r, b)                % Graspling fixture conic constraint
         ];                                                           

    % Equality constraints
%     omega(:,end) = QuaternionAlgebra.RotateVector(QuaternionAlgebra.quaternion_inverse( beta(1:4,1) ), omega(:,end));
%     omega_tf =     QuaternionAlgebra.RotateVector(QuaternionAlgebra.quaternion_inverse( q_t(:,end) ),  omega_t(:,end));

    ceq = [
            dot(s(1:4,:), s(1:4,:), 1).'-1; ...         % Quaternion norm constraint
%             omega(:,end) - omega_tf ...                 % Relative angular velocity
          ]; 
end