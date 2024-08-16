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
    cos_l = cos(t); 
    sin_l = sin(t);
    k = 2 ./ ( 1 - dot(C(4:5,:), C(4:5,:), 1) ) .* ( C(5,:) .* cos_l - C(4,:) .* sin_l );
    r_h = delta ./ w;

    % Pre-allocation 
    B = zeros(6, 3 * size(C,2));

    % Compute the control input matrix
    for i = 1:size(C,2)
        % Semilatus rectum / angular momentum
        B(1, 3 * (i-1) + 2) = 2 * r_h(i) * C(1,i) ;

        % Projection of the eccentricities
        B(2, 3 * (i-1) + 1) = +delta(i) * sin_l(i);
        B(3, 3 * (i-1) + 1) = -delta(i) * cos_l(i);

        B(2, 3 * (i-1) + 2) = r_h(i) * ( (w(i) + 1) * cos_l(i) + C(2,i) );
        B(3, 3 * (i-1) + 2) = r_h(i) * ( (w(i) + 1) * sin_l(i) + C(3,i) );

        B(2, 3 * (i-1) + 3) = +r_h(i) * k(i) * C(3,i);
        B(3, 3 * (i-1) + 3) = -r_h(i) * k(i) * C(2,i);

        % MRPs
        omega = r_h(i) * [cos_l(i); sin_l(i); k];
        omega = 0.25 * QuaternionAlgebra.Quat2Matrix([C(4:5,i); 0; -1]) * omega;
        B(4:5, 3 * (i-1) + 3) = omega(1:2,1);

        % Longitude / ideal latitude
        B(6, 3 * (i-1) + 3) = r_h(i) * (-k(i));
    end
end