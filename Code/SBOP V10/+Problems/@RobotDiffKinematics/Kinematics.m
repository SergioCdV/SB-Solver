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

function [T, J] = Kinematics(n, type, transformation, s)    
    % Preallocation 
    P = zeros(4, n);        % Third column of the homogeneous Euclidean transformation
    Z = zeros(3, n);        % Third column of the rotation matrix of step n
    J = zeros(6, n);        % Vector of n vectors of dimension 6 x 1 
    
    % Compute the z and p vectors  
    Z(:,1) = [0;0;1];
    P(:,1) = [0;0;0;1];

    for i = 1:n
        % Compute the transformation matrix of the n joint 
        [A, R] = transformation(i, s);

        % Perform the transformation 
        Z(:,i+1) = R * Z(:,i);
        P(:,i+1) = A * P(:,i);
    end

    T = A;

    % Assemble the Jacobian
    for i = 1:n
        % Compute the state variable 
        z = Z(:,i);
        pe = P(1:3,end);
        pm = P(1:3,i);

        % Assemble the column vector
        if (type == 1)
            J(:,i) = [cross(z, pe-pm); z];
        else
            J(:,i) = [z; zeros(3,1)];
        end
    end
end