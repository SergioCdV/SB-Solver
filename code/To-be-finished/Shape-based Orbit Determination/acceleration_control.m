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
    % Compute the control vector
    r = sqrt(C(1,:).^2+C(3,:).^2);

    % Normalizing factor
    switch (method)
        case 'Sundman'
            h = angular_momentum(C); 
            eta = r.^2./h;
            c = eta;
        otherwise
            c = tf;
    end

    % Compute the control vector
    u = [C(7,:)-C(1,:).*C(5,:).^2+c.^2.*mu.*C(1,:)/r.^3; ...
         C(1,:).*C(8,:)+2*C(4,:).*C(5,:); ... 
         C(9,:)+c.^2.*mu.*C(3,:)/r.^3];
end
