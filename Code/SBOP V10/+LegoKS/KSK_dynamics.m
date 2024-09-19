
function [ds] = KSK_dynamics(E, t, s)
    % Factorization
    s = reshape(s, 8, []);      % Each point in the fiber

    % State variables 
     u = s(1:4,:);              % State on the fiber
    du = s(5:8,:);              % Derivative of the fiber

    % Frequency of the oscillation
    omega = sqrt( -E/2 );   

    % Dynamics 
    ds = zeros(size(s));
    ds(1:4,:) = du;
    ds(5:8,:) = -omega.^2.' .* u;
end