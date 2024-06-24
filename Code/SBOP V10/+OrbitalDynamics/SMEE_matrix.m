%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 22/11/2022

%% Control input matrix %% 
% Function to transform classical orbital elements into Cartesian or
% viceversa

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the orbital equinoctial state vector to be transformed 

% Outputs: - matrix B, the control input matrix

function [B] = SMEE_matrix(mu, t, C)
    % Compute the auxiliary variables
    w = 1 + C(2,:).* cos(t) + C(3,:) .* sin(t);
    delta = sqrt(C(1,:) ./ mu);

    % Pre-allocation 
    B = zeros(6, 3 * size(C,2));

    % Compute the control input matrix
    for i = 1:size(C,2)
        % Constants 
        k = 2 / ( 1 - dot(C(4:5,i), C(4:5,i), 1) ) * ( C(5,i) .* cos(t(i)) - C(4,i) .* sin(t(i)) );

        % Semilatus rectum / angular momentum
        B(1, 3 * (i-1) + 2) = 2 * delta(i) * C(1,i) / w(i);

        % Projection of the eccentricities
        B(2, 3 * (i-1) + 1) = +delta(i) * sin(t(i));
        B(3, 3 * (i-1) + 1) = -delta(i) * cos(t(i));

        B(2, 3 * (i-1) + 2) = delta(i) / w(i) * ( (w(i) + 1) * cos(t(i)) + C(2,i) );
        B(3, 3 * (i-1) + 2) = delta(i) / w(i) * ( (w(i) + 1) * sin(t(i)) + C(3,i) );

        B(2, 3 * (i-1) + 3) = -delta(i) * C(3,i) / w(i) * (-k);
        B(3, 3 * (i-1) + 3) = +delta(i) * C(2,i) / w(i) * (-k);

        % MRPs
        omega = delta(i) / ( 2 * w(i) ) * [cos(t(i)); sin(t(i)); k];
        omega = 0.25 * QuaternionAlgebra.Quat2Matrix([C(4:5,i); 0; -1]) * omega;
        B(4:5, 3 * (i-1) + 3) = omega(1:2,1);

        % Longitude / ideal latitude
        B(6, 3 * (i-1) + 3) = -delta(i) * k / w(i);
    end
end