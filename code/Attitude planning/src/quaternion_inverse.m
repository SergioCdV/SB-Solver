%% Project: 
% Date: 16/04/22

%% Quaternion inverse %%
% This file contains the function to compute the inverse of a given quaternion

% Inputs: - quaternion q1, to be inverted. 

% Ouputs: - quaternion q, the output of the quaternion inversion.

% New version updates: 

function [q] = quaternion_inverse(q1)
    % Initialize the quaternion 
    q = zeros(4,1);           

    % Perform the quaternion conjugation 
    q(2:4) = -q1(2:4);
    q(1) = q1(1);     
    
    % Compute the quaternion inverse 
    q = q/norm(q1)^2;
end