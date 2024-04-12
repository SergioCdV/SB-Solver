

function [A] = Quat2Matrix(q)
    % Correct dimensions 
    q = q.';

    % Rotation matrices
    A = (q(:,4).^2 - dot(q(:,1:3), q(:,1:3), 2)) * repmat(eye(3), size(q,1), 1) + 2 * q(:,1:3).' * q(:,1:3) - 2 * q(:,4) * QuaternionAlgebra.hat_map(q(:,1:3).');
end