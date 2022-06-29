%% Project: 
% Date: 16/04/22

%% Quaternion product %%
% This file contains the function to compute the quaternion product between two quaternions following Hamilton's
% convention

% Inputs: - quaternion q1. 
%         - quaternion q2. 

% Ouputs: - quaternion q, the output of the quaternion products.

% New version updates: 

function [q] = quaternion_product(q1, q2)
    %Initialize the quaternion 
    q = zeros(4,1);                  
    
    %Compute the quaternion product 
    q(2:4) = q1(1)*q2(2:4)+q2(1)*q1(2:4)+cross(q1(2:4),q2(2:4)); 
    q(1) = q1(1)*q2(1)-dot(q1(2:4),q2(2:4));
end