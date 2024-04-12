%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Inequality constraints
%     detJ = zeros(1,size(tau,2));
%     objects = [];
% 
%     if (~isempty(objects))
%         dm = zeros(1, ( factorial(length(objects)) / (factorial(length(objects)-2) * 2)) * size(tau,2));
%     else
%         dm = [];
%     end
% 
%     counter = 1;
%     v = zeros(6,length(tau));
%     z_pos = zeros(1,length(tau));
%     for i = 1:size(tau,2)
%         [T, J] = Problems.RobotDiffKinematics.Kinematics(obj.StateDim, ...
%                                                  @(j,s)Problems.RobotDiffKinematics.ur3_dkinematics(obj, j, s), ...
%                                                  s(:,i));
%         % Jacobian singularities
%         detJ(i) = det(J)^2;
% 
%         v(:,i) = J * s(7:12,i);
% 
%         z_pos(i) = T(3,4);
% 
%         % Compute the collisions
%         for j = 1:length(objects)
%             for k = j+1:length(objects)
%                 index_1 = objects(j);
%                 index_2 = objects(k);
% 
%                 r1 = T(1:3, 4 * index_1);
%                 A1 = T(1:3, 1 + 4 * (index_1-1): 3 + 4 * (index_1-1));
%                 r2 = T(1:3, 4 * index_2);
%                 A2 = T(1:3, 1 + 4 * (index_2-1): 3 + 4 * (index_2-1));
% 
%                 v1 = obj.body_vertex(index_1);
%                 v2 = obj.body_vertex(index_2);
% 
%                 dm(counter) = obj.collision_constraint(r1, A1, v1, r2, A2, v2, 1);
%                 counter = counter + 1;
%             end
%         end
%     end
% 
%     epsilon = params(7);
% 
    DH_parameters.base =   reshape(params(8:10), 1, 3).'; 
    DH_parameters.theta =  reshape(params(11:16), 1, obj.StateDim).';
    DH_parameters.alpha =  reshape(params(17:22), 1, obj.StateDim).';
    DH_parameters.offset = reshape(params(23:28), 1, obj.StateDim).';
    DH_parameters.a =      reshape(params(29:34), 1, obj.StateDim).';
    DH_parameters.d =      reshape(params(35:40), 1, obj.StateDim).';
    DH_parameters.type =   reshape(params(41:46), 1, obj.StateDim).';

    c = [
         -reshape(s(7:9,:), 1, [])-params(3) ...        % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
         +reshape(s(7:9,:), 1, [])-params(3) ...        % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
         -reshape(s(10:12,:), 1, [])-params(4) ...      % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
         +reshape(s(10:12,:), 1, [])-params(4) ...      % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
%          -reshape(v(1:3,:), 1, [])-params(5) ...        % Constraint on the frame linear velocity (infinity norm in epigraph form)
%          +reshape(v(1:3,:), 1, [])-params(5) ...        % Constraint on the frame linear velocity (infinity norm in epigraph form)
%          -reshape(v(4:6,:), 1, [])-params(6) ...        % Constraint on the frame angular velocity (infinity norm in epigraph form)
%          +reshape(v(4:6,:), 1, [])-params(6) ...        % Constraint on the frame angular velocity (infinity norm in epigraph form)
%          +dm ...                                        % Collisions constraints
%          epsilon-reshape(detJ, 1, []) ...               % Singularity constraint
%          -z_pos                                         % Constraint to lie above the table
%          reshape(s(1:6,1)-s(1:6,end), 1, [])-pi ...     % Constraint on the frame angular velocity (infinity norm in epigraph form)
%          reshape(s(1:6,end)-s(1:6,1), 1, [])-pi ...     % Constraint on the frame angular velocity (infinity norm in epigraph form)
        ];   
% 
%     % Equality constraints    
%     if (length(params) < 21)
%         s_ref = reshape(params(8:20), [], 1);           % Desired waypoint
% 
%         q_f(4,1) = sqrt(1 + T(1,1) + T(2,2) + T(3,3))/2;
%         q_f(1,1) = (T(3,2) - T(2,3)) / (4 * q_f(4));
%         q_f(2,1) = (T(1,3) - T(3,1)) / (4 * q_f(4));
%         q_f(3,1) = (T(2,1) - T(1,2)) / (4 * q_f(4));
%     
%         dq = QuaternionAlgebra.right_isoclinic( q_f ) * QuaternionAlgebra.quaternion_inverse( s_ref(4:7,1) );
% 
%         ceq = [
%                T(1:3,end)-s_ref(1:3,1); ...            % Final linear position constraint
%                J*s(7:12,end)-s_ref(8:13,1); ...        % Final linear and angular velocity constraint
%                dq(4) - 1;                              % Final attitude constraint
%                ];   
%     else
%         ceq = [];
%     end

      ceq = [];
end