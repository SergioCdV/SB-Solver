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
    f = -[-c.*C(4:6,:); c.*mu.*C(1,:)./r.^3; zeros(1,size(C,2)); c.*mu.*C(3,:)./r.^3];     % Dynamics vector in first order form

    % Compute the acceleration term
    a = zeros(6,length(t));
    for i = 1:length(t)
        v = sum(D(i,:)/tf.*C(1:3,:),2);
        gamma = sum(D(i,:)/tf.*C(4:6,:),2);
        a(:,i) = [v; gamma(1,1)-C(1,1).*v(2,1).^2; C(1,1).*gamma(2,1)+2*v(1,1).*v(2,1); gamma(3,1)];
    end

    % Compute the dynamics residual
    res = f+[zeros(3,length(t));u]-a;
end
