%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 22/11/2022

%% Control input matrix %% 
% Function to transform classical orbital elements into Cartesian or
% viceversa

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the orbital equinoctial state vector to be transformed 

% Outputs: - matrix B, the control input matrix

function [B] = MEE_matrix(mu, t, C)
    % Compute the auxiliary variables
    w = 1+C(2,:).*cos(t)+C(3,:).*sin(t);
    s = 1+C(4,:).^2+C(5,:).^2;

    % Compute the control input matrix
    for i = 1:size(C,2)
        B(1,2) = 2*C(1,i)./w(i).*sqrt(C(1,i)./mu);
        B(2,1) = sin(t).*sqrt(C(1,i)./mu);
        B(2,2) = sqrt(C(1,i)/mu)./w(i).*((w(i)+1).*cos(t)+C(2,i));
        B(2,3) = -sqrt(C(1,i)/mu).*C(3,i)./w(i).*(C(4,i).*sin(t)-C(5,i).*cos(t));
        B(3,1) = -cos(t).*sqrt(C(1,i)/mu);
        B(3,2) = sqrt(C(1,i)/mu)./w(i).*((w(i)+1).*sin(t)+C(3,i));
        B(3,3) = sqrt(C(1,i)/mu).*C(2,i)./w(i).*(C(4,i).*sin(t)-C(5,i).*cos(t));
        B(4,3) = sqrt(C(1,i)/mu).*s(i).*cos(t)./(2*w(i));
        B(5,3) = sqrt(C(1,i)/mu).*s(i).*sin(t)./(2*w(i));
        B(6,3) = sqrt(C(1,i)/mu)./w(i).*(C(4,i).*sin(t)-C(5,i).*cos(t));
    end
end