%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    R = dot(params(4:6).'-s(1:3,end), params(4:6).'-s(1:3,end), 1)-params(3)^2; 

    c = [dot(u(1:3,:),u(1:3,:),1)-params(2)^2 ...  % Constraint on the force magnitude (second order cone)
         dot(u(4:6,:),u(4:6,:),1)-params(22)^2 ... % Constraint on the torque magnitude (second order cone)
         dot(s(4:7,:),s(4:7,:),1)-1 ...            % Quaternion norm constraint
         R].';                                     % Graspling inequality

    % Equality constraints
    ceq = s(8:10,end)-cross(params(7:9).', s(1:3,end)-params(10:12).');
    ceq = reshape(ceq, [], 1);
end