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
    c = tf;                                                                                % Normalizing factor
    f = [C(4:6,:); -mu.*C(1,:)./r.^3; zeros(1,size(C,2)); -mu.*C(3,:)./r.^3];    % Dynamics vector in first order form

    % Compute the acceleration term
    a = zeros(6,length(t));
    v = C(1:3,:)*D.';
    gamma = C(4:6,:)*D.';
    for i = 1:length(t)
        a(:,i) = [v(:,i); gamma(1,i)-C(1,i).*v(2,i).^2; C(1,i).*gamma(2,i)+2*v(1,i).*v(2,i); gamma(3,i)];
    end

    % Compute the dynamics residual
    res = a-tf*(f+[zeros(3,length(t));u]);
end
