%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the 9xm state vector 
%         - scalar tf, the final time of flight
%         - vector t, the independent variable evolution

% Outputs: - vector u, the nondimensional 3xm control vector
%          - vector dv, the inertial velocity field 
%          - vector f, the dynamics vector field

function [u, dv, f] = acceleration_control(mu, C, tf, t)
    % Location of the primaries 
    R(:,1) = [-mu; 0; 0];           % Location of the first primary
    R(:,2) = [1-mu; 0; 0];          % Location of the second primary

    r(1:3,:) = C(1:3,:)-R(:,1);     % Relative position to the first primary
    r(4:6,:) = C(4:6,:)-R(:,2);     % Relative position to the second primary

    % Linear terms of the equations of motion
    c = tf;                                                                                                          % Normalizing factor
    f = -(1-mu)*r(1:3,:)./sqrt(dot(r(1:3,:),r(1:3,:),1)).^3-(mu)*r(4:6,:)./sqrt(dot(r(4:6,:),r(4:6,:),1)).^3;        % Acceleration vector
    a = C(7:9,:);                                                                                                    % Inertial acceleration field
    dv = C(4:6,:);                                                                                                   % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = a-c^2*f;
end
