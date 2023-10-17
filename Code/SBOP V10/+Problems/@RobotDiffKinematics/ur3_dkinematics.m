function [T, R] = ur3_dkinematics(obj, base, theta, alpha, offset, a, d, type, n, q)
    % Preallocation
    T = eye(4);

    for i = 1:n
        % Assemble the state vector
        if (type(i) == 1)
            theta(i) = q(i);
        elseif (type(i) == 0)
            d(i) = q(i);
        end

        % Compute the homogeneous matrix
        ct = cos(theta(i) + offset(i));
        st = sin(theta(i) + offset(i));
        ca = cos(alpha(i));
        sa = sin(alpha(i));
        
        T = [ct -st*ca st*sa a(i)*ct ; ...
             st ct*ca -ct*sa a(i)*st ; ...
             0 sa ca d(i); ...
             0 0 0 1] * T;
        T(1:3,4) = T(1:3,4) + base;
    end
    
    % Extract the rotation matrix
    R = T(1:3,1:3);
end