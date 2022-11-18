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
    f = zeros(N,10);                                                     % Vector field

    % Dynamics
    f(:,1:5) = s(:,6:10)/tf;                                             % Complete dynamics
    for i = 1:N
        r = dot(s(i,1:4), s(i,1:4));
        L = KS_matrix(s(i,1:4)); 
        f(i,6:10) = [-mu*s(i,5)/4.*s(i,1:4)+r*u(1:4,i).'*L/2 -(4/mu)*dot(s(i,6:9)/tf,u(1:4,i).'*L)];           
    end

    ds = tf*f;
end