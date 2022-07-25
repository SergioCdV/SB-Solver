%% Project: 
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the 9xm state vector 
%         - scalar tf, the final time of flight
%         - string method, the collocation distribution to be used

% Outputs: - vector u, the nondimensional 3xm control vector

function [u] = acceleration_control(mu,C,tf,method)
    % Compute the radius vector
    r = sqrt(C(1,:).^2+C(3,:).^2);

    % Compute the control vector as a residual of the dynamics
    switch (method)
        case 'Regularized'
            % Normalizing factor
            c = tf;

            % Derivative of the radius with the arclength 
            dr = dot(C(1,:).*C(4,:)+C(3,:).*C(6,:))./r;

            % Compute the control vector as a dynamics residual
            u = [C(7,:)-dr.*C(4,:)./r+c.^2.*mu.*C(1,:)./r-C(1,:).*C(5,:).^2; ...
                 C(1,:).*C(8,:)+2*C(4,:).*C(5,:); ... 
                 C(9,:)-dr.*C(6,:)./r+c.^2.*mu.*C(3,:)./r];

            u(2,:) = u(2,:)./r.^2; 

        otherwise
            % Normalizing factor
            c = tf;

            % Compute the control vector as a dynamics residual
            u = [C(7,:)-C(1,:).*C(5,:).^2+c.^2.*mu.*C(1,:)./r.^3; ...
                 C(1,:).*C(8,:)+2*C(4,:).*C(5,:); ... 
                 C(9,:)+c.^2.*mu.*C(3,:)./r.^3];
    end
end
