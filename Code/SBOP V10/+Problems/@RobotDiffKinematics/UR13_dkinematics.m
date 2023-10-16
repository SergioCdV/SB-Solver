%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/10/2023

%% UR13 direct kinematics %% 
% Function to compute the transformation matrix of a robot from the n-1 system to the n one

% Inputs: - scalar n, the joint index
%         - vector s, the current joint state vector, of dimensions n x 1

% Outputs: - matrix T, the end-effector to base homogeneous transformation matrix

function [T] = UR13_dkinematics(n, s)   
    % Compute the transformation matrices
    switch (n)
        case 1
        case 2 
        case 3 
        case 4 
        case 5 
        case 6
        otherwise
            error('An invalid joint index was selected...');
    end
end