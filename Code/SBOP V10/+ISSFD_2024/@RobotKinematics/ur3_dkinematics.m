function [T] = ur3_dkinematics(DH_parameters, i, q)
    % Constants 
    base =   DH_parameters.base;
    theta =  DH_parameters.theta;
    alpha =  DH_parameters.alpha;
    offset = DH_parameters.offset;
    a =      DH_parameters.a;
    d =      DH_parameters.d;
    type =   DH_parameters.type;

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