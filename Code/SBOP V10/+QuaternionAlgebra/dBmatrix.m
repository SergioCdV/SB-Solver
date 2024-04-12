
%% Auxiliary function 
function [Bdot] = dBmatrix(q, dq)
    Bdot = zeros(3,3);
    s = -2 * dot(q, dq);
    Bdot(1,1) = s + 4 * (q(1) * dq(1));
    Bdot(1,2) = 2 * (-dq(3) + q(1) * dq(2) + dq(1) * q(2));
    Bdot(1,3) = 2 * ( dq(2) + q(1) * dq(3) + dq(1) * q(3));
    Bdot(2,1) = 2 * ( dq(3) + q(1) * dq(2) + dq(1) * q(2));
    Bdot(2,2) = s + 4 * (q(2) * dq(2));
    Bdot(2,3) = 2 * (-dq(1) + q(2) * dq(3) + dq(2) * q(3));
    Bdot(3,1) = 2 * (-dq(2) + q(1) * dq(3) + dq(1) * q(3));
    Bdot(3,2) = 2 * ( dq(1) + q(2) * dq(3) + dq(2) * q(3));
    Bdot(3,3) = s + 4 * (q(3) * dq(3));
end