%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - scalar mu, the gravitational parameter of the system
%         - matrix D, the differentiation matrix
%         - array C, the 6xm state vector 
%         - array u, the 3xm control vector 
%         - scalar tf, the final time of flight
%         - vector t, the independent variable evolution

% Outputs: - vector res, the dynamics residual

function [res] = acceleration_control(D, C, u, tf, t)
    % Linear terms of the equations of motion
    f = [C(2,:); zeros(1,length(t))];    % Dynamics vector in first order form

    % Compute the acceleration term
    dC = C*D.';
    a = dC;

    % Compute the dynamics residual
    res = a-tf/2*(f+[zeros(1,length(t)); u]);
end
