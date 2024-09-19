
function [R] = GroupAction(beta)
    % Pre-computation
    sin_beta = sin(beta); 
    cos_beta = cos(beta);

    % Rotation matrix
    R = [cos_beta 0 0 -sin_beta; ... 
        0 cos_beta sin_beta 0; ...
        0 -sin_beta cos_beta 0; ...
        sin_beta 0 0 cos_beta];
end