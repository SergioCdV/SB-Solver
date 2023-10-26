%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    detJ = zeros(1,size(tau,2));

    for i = 2:size(tau,2)
        [T, J] = Problems.RobotDiffKinematics.Kinematics(obj.StateDim, ...
                                                 @(j,s)Problems.RobotDiffKinematics.ur3_dkinematics(obj, j, s), ...
                                                 s(:,i));
        detJ(i) = det(J)^2;
    end

    c = [-reshape(s(7:9,:), 1, [])-params(3) ...   % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
          reshape(s(7:9,:), 1, [])-params(3) ...   % Constraint on the angular velocity           (infinity norm in epigraph form)
         -reshape(s(10:12,:), 1, [])-params(4) ... % Constraint on the angular acceleration       (infinity norm in epigraph form)
          reshape(s(10:12,:), 1, [])-params(4) ... % Constraint on the angular acceleration       (infinity norm in epigraph form)
          % 1e-6-reshape(detJ, 1, []);
        ];                                                                                      

    % Equality constraints
    s_ref = reshape(params(5:end), [], 1);
    ceq = [% T(1:3,end)-s_ref(1:3,1); ...
           % reshape(T(1:3,end-3:end-1)-QuaternionAlgebra.Quat2Matrix(s_ref(4:7,1)), [], 1); ...
           % J*s(7:12,end)-s_ref(8:13,1);
           ];    
end