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

function [res] = acceleration_control(mu, D, C, u, tf, t)
    % Compute the radius vector
    r = sqrt(C(1,:).^2+C(3,:).^2);

    % Linear terms of the equations of motion
    f = [C(4:6,:); -mu.*C(1,:)./r.^3; zeros(1,size(C,2)); -mu.*C(3,:)./r.^3];    % Dynamics vector in first order form

    % Compute the acceleration term
    dC = C*D.';
    a = [dC(1:3,:); dC(4,:)-C(1,:).*C(5,:).^2; C(1,:).*dC(5,:)+2*C(4,:).*C(5,:); dC(6,:)];

    % Compute the dynamics residual
    res = a-tf*(f+[zeros(3,length(t)); u]);
end
