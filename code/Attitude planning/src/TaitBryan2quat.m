%% Project: Shape-based attitude planning %%
% Date: 29/06/22

%% Quaternion to Tait-Bryan angles %%
% This file contains the function to generate a rotation matrix from a
% general quaternion to the instantenous values of the associated Tait-Bryan angles. 

% Inputs: - argument input, the general quaternion or the associated Tait-Bryan angles
%         - boolean quaternion_flag, to indicate whether the input is a
%           quaternion or a rotation matrix

% Ouputs: - vector output, containing the roll, pitch and yaw angles or an
%           associated quaternion together with the associated angular velocity

% New version updates: 

function [output] = TaitBryan2quat(quaternion_flag, input)
    % Branch the algorithm 
    if (quaternion_flag)
        % Associated Tait-Bryan angles
        q = input;                                                      % Associated quaternion
        yaw = atan2(q(4)*q(3)+q(1)*q(2), 0.5-(q(2)^2+q(3)^2));          % Yaw angle
        pitch = asin(2*(q(4)*q(2)-q(3)*q(1)));                          % Pitch angle
        roll = atan2(q(4)*q(1)+q(2)*q(3), 0.5-(q(1)^2+q(2)^2));         % Roll angle
        
        % Final Tait-Bryan angles
        output = [yaw pitch roll input(5:8)];  

    else
        % Associated quaternion
        yaw = input(1);               % Yaw angle 
        pitch = input(2);             % Pitch angle 
        roll = input(3);              % Roll angle
        
        q(1) = sin(roll/2)*cos(pitch/2)*cos(yaw/2)-cos(roll/2)*sin(pitch/2)*sin(yaw/2);   
        q(2) = cos(roll/2)*sin(pitch/2)*cos(yaw/2)+sin(roll/2)*cos(pitch/2)*sin(yaw/2); 
        q(3) = cos(roll/2)*cos(pitch/2)*sin(yaw/2)-sin(roll/2)*sin(pitch/2)*cos(yaw/2); 
        q(4) = cos(roll/2)*cos(pitch/2)*cos(yaw/2)+sin(roll/2)*sin(pitch/2)*sin(yaw/2);  
  
        % Final quaternion 
        output = [q input(4:6)];   
    end
end