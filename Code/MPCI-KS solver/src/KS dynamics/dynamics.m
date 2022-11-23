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
    N = size(u,2);                                                      % Number of evaluation points
    f = zeros(N,8);                                                     % Vector field

    % Dynamics
    f(:,1:4) = s(:,5:8);                                             % Complete dynamics
    r = dot(s(:,1:4), s(:,1:4), 2);
    E = -2*(mu/2-dot(s(:,5:8),s(:,5:8),2))./r;
    for i = 1:N
        L = KS_matrix(s(i,1:4)); 
        f(i,5:8) = -E(i)/2.*s(i,1:4)+tf*r(i)*u(1:4,i).'*L/2;           
    end

    ds = f;
end