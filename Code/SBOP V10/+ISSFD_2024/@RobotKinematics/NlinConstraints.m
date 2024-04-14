%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, tau, s, u)
    % Constants 
    % DH parameters 
    DH_parameters.base =   reshape(params(8:10), 1, 3).'; 
    DH_parameters.theta =  reshape(params(11:16), 1, obj.StateDim).';
    DH_parameters.alpha =  reshape(params(17:22), 1, obj.StateDim).';
    DH_parameters.offset = reshape(params(23:28), 1, obj.StateDim).';
    DH_parameters.a =      reshape(params(29:34), 1, obj.StateDim).';
    DH_parameters.d =      reshape(params(35:40), 1, obj.StateDim).';
    DH_parameters.type =   reshape(params(41:46), 1, obj.StateDim).';

    % Inequality constraints
%     objects = [];
% 
%     if (~isempty(objects))
%         dm = zeros(1, ( factorial(length(objects)) / (factorial(length(objects)-2) * 2)) * size(tau,2));
%     else
%         dm = [];
%     end
% 
%     counter = 1;
    detJ = zeros(1, size(tau,2));
    v = zeros(6, size(tau,2));

%     z_pos = zeros(1,length(tau));
    for i = size(tau,2):size(tau,2)
        [T, J] = ISSFD_2024.RobotKinematics.Kinematics(obj.StateDim, @(j,s)ISSFD_2024.RobotKinematics.ur3_dkinematics(DH_parameters, j, s), s(:,i));
        
        % Jacobian singularities
        detJ(i) = det(J)^2;
        
        % Velocity of the end frame
        v(:,i) = J * s(7:12,i);

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
    end
            
    % Reference position and velocity 
    R = T(1:3,end-2:end);                   % Final rotation matrix
    q = QuaternionAlgebra.Matrix2Quat(R);   % Associated quaternion

    sigma = QuaternionAlgebra.MPR2Quat(1, 1, q, false);

    r_ref = params(47:49).';                % Reference position
    sigma_ref = params(50:52).';            % Reference attitude
    v_ref = params(53:58).';                % Reference velocity

    r_eff = T(1:3,end);                     % End effector position 
    res(1:3,1) = r_eff - r_ref;             % Error to the reference position
    res(4:6,1) = sigma - sigma_ref;         % Error to the reference attitude 
    res(7:12,1) = v(:,end) - v_ref;         % Error to the velocities

    % Inequalities
    c = [
         +reshape(u(1:4,:), 1, [])-params(3) ...        % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
         -reshape(u(1:4,:), 1, [])-params(3)  ...       % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
         +reshape(u(5:6,:), 1, [])-params(4) ...        % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
         -reshape(u(5:6,:), 1, [])-params(4) ...        % Constraint on the angular velocity magnitude (infinity norm in epigraph form)
         reshape(+res -1E-3 , 1, []) ...                % Relaxed equalities
         reshape(-res -1E-3 , 1, []) ...                % Relaxed equalities
%          reshape(+v(1:3,:)-params(5), 1, []) ...        % Constraint on the frame linear velocity (infinity norm in epigraph form)
%          reshape(-v(1:3,:)-params(5), 1, []) ...        % Constraint on the frame linear velocity (infinity norm in epigraph form)
%          reshape(+v(4:6,:)-params(6), 1, []) ...        % Constraint on the frame angular velocity (infinity norm in epigraph form)
%          reshape(-v(4:6,:)-params(6), 1, []) ...        % Constraint on the frame angular velocity (infinity norm in epigraph form)
%          params(7)-reshape(detJ, 1, []) ...             % Singularity constraint
        ];   

    ceq = [];
end