%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/11/22

%% Dynamics %%
% Function to compute the dynamics vector field to be integrated in the
% MCPI method

% Inputs: - scalar mu, the gravitational parameter of the system
%         - scalar tf, the final time of flight
%         - array u, the mx3 control vector 
%         - vector tau, the time domain
%         - array s, the mx6 state vector 

% Outputs: - vector ds, the dynamics vector field


% State dynamics 
function [ds] = dynamics(mu, tf, u, tau, s)
    N = size(u,2);                                                       % Number of evaluation points

    s(:,4:6) = s(:,4:6)/tf;
    
    r = sqrt(dot(s(:,[1 3]), s(:,[1 3]), 2));                            % Radial distance
    a = [s(:,1).*s(:,5).^2 -2*s(:,4).*s(:,5)./s(:,1) zeros(N,1)];        % Inertial acceleration terms
    f = -[mu.*s(:,1)./r.^3 zeros(N,1) mu.*s(:,3)./r.^3];                 % Gravity field
    f = [s(:,4:6) a+f];                                                  % Complete dynamics
    ds = tf*(f+[zeros(size(u,2),3) u.']);                                % Controlled plant
end