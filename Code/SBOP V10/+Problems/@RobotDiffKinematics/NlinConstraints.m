%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
    detJ = zeros(1,size(tau,2));
    objects = [1 2 6];

    if (~isempty(objects))
        dm = zeros(1, factorial(length(objects)) / (factorial(length(objects)-2) * 2) * size(tau,2));
    end

    counter = 1;
    for i = 1:size(tau,2)
        [T, J] = Problems.RobotDiffKinematics.Kinematics(obj.StateDim, ...
                                                 @(j,s)Problems.RobotDiffKinematics.ur3_dkinematics(obj, j, s), ...
                                                 s(:,i));

        % Jacobian singularities
        detJ(i) = det(J)^2;

        % Compute the collisions
        for j = 1:length(objects)
            for k = j+1:length(objects)

                index_1 = objects(j);
                index_2 = objects(k);

                r1 = T(1:3, 4 * index_1);
                A1 = T(1:3, 1 + 4 * (index_1-1): 3 + 4 * (index_1-1));
                r2 = T(1:3, 4 * index_2);
                A2 = T(1:3, 1 + 4 * (index_2-1): 3 + 4 * (index_2-1));

                v1 = rand(3,4);
                v2 = rand(3,4);

                dm(counter) = obj.collision_constraint(r1, A1, v1, r2, A2, v2, 10);
                counter = counter + 1;
            end
        end
    end

    epsilon = params(5);

    c = [
         -reshape(s(7:9,:), 1, [])-params(3) ...        % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
         +reshape(s(7:9,:), 1, [])-params(3) ...        % Constraint on the angular velocity           (infinity norm in epigraph form)
         -reshape(s(10:12,:), 1, [])-params(4) ...      % Constraint on the angular acceleration       (infinity norm in epigraph form)
         +reshape(s(10:12,:), 1, [])-params(4) ...      % Constraint on the angular acceleration       (infinity norm in epigraph form)
          epsilon-reshape(detJ, 1, []) ...              % Singularity constraint 
%           dm/abs(max(dm)) ...                                        % Collisions constraints
        ];   

    max(c)

    % Equality constraints
    s_ref = reshape(params(6:end), [], 1);                                                      % Desired waypoint
    ceq = [
%            T(1:3,end)-s_ref(1:3,1); ...                                                         % Final linear position constraint
           J*s(7:12,end)-s_ref(8:13,1); ...                                                     % Final linear and angular velocity constraint
%            reshape(T(1:3,end-3:end-1)-QuaternionAlgebra.Quat2Matrix(s_ref(4:7,1)), [], 1); ...  % Final attitude constraint
           ];    
end