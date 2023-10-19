%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/10/2023

%% Kinematics %% 
% Function to compute the Jacobian and the transformation matrix of a robot
% from the base to the end-effector 

% Inputs: - scalar n, the joint space dimension
%         - vector type, indicating the type of the nth joint (either revolute (0) or prismatic (1))
%         - function handle transformation, containing the different homogeneous transformation matrices for each joint
%         - vector s, the current joint state vector, of dimensions n x 1

% Outputs: - matrix T, the end-effector to base homogeneous transformation matrix
%          - matrix J, the differential kinematics Jacobian 

function [T, J] = Kinematics(n, transformation, s)    
    % Preallocation 
    d = zeros(3, n+1);      % Third column of the homogeneous Euclidean transformation
    Z = zeros(3, n);        % Third column of the rotation matrix of step n
    J = zeros(6, n);        % Vector of n vectors of dimension 6 x 1 
    
    % Compute the z and p vectors  
    Z(:,1) = [0;0;1];

    T = zeros(4, (n+1)*4);
    T(1:4,1:4) = eye(4);

    for i = 1:n
        % Compute the transformation matrix of the n joint 
        [A] = transformation(i, s(i));
        T(:, 1+4*i:4*(i+1)) = T(:, 1+4*(i-1):4*i) * A;
    end

    for i = 1:n
        % Perform the transformation 
        Z(:,i) = T(1:3,1+4*(i-1):3+4*(i-1)) * [0; 0; 1];
        d(:,i+1) = T(1:3,4*i);
    end

    % Assemble the Jacobian
    for i = 1:n
        % Compute the state variable 
        z = Z(:,i);
        D = d(:,n+1)-d(:,i);

        % Assemble the column vector
        if (1)
            J(:,i) = [cross(z, D); z];
        else
            J(:,i) = [z; zeros(3,1)];
        end
    end
end