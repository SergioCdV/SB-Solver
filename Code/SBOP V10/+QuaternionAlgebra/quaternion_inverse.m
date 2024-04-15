

function [qi] = quaternion_inverse(q)
    qi = q;
    qi(1:3,:) = -qi(1:3,:) ./ sqrt(dot(qi, qi, 1));
end