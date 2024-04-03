
function [Q] = ECI2LVLH(r, v, direction)

    % Compute the rotation matrix 
    h = cross(r,v);                             % Angular momentum vector

    o21 = -h ./ sqrt(dot(h,h,1));               % LVLH y+ definition
    o31 = -r(1:3,:) ./ sqrt(dot(r,r,1));        % LVLH z+ definition
    o11 = cross(o21, o31);                      % LVLH x+ definition

    Q = reshape([o11; o21; o31], 3, []).';      % Rotation matrix from the LVLH to the ECI frames

    if (~direction)
        Q = Q.';
    end

end