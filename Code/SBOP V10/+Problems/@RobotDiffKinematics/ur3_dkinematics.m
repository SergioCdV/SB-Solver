function [T] = ur3_dkinematics(obj, i, q)
    % Constants 
    base = [0; 0; 0];
    theta = [0 0 0 0 0 0];
    alpha = [pi/2 0 0 pi/2 -pi/2 0];
    offset = [0 0 0 0 0 0];
    a = [0 -0.24365 -0.21325 0 0 0];
    d = [0.15185 0 0 0.1124 0.08535 0.0921];
    type = ones(6,1);

    % Assemble the state vector
    if (type(i) == 1)
        theta = q;
    elseif (type(i) == 0)
        d = q;
    end

    % Compute the homogeneous matrix
    ct = cos(theta + offset(i));
    st = sin(theta + offset(i));
    ca = cos(alpha(i));
    sa = sin(alpha(i));
    
    T = [ct -st*ca +st*sa a(i)*ct; ...
         st +ct*ca -ct*sa a(i)*st; ...
          0     sa     ca d(i); ...
          0      0      0 1];
    T(1:3,4) = T(1:3,4) + base;
end