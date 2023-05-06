%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 22/11/2022

%% Control input matrix %% 
% Function to transform classical orbital elements into Cartesian or
% viceversa

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array s, the orbital equinoctial state vector to be transformed 

% Outputs: - matrix B, the control input matrix

function [B] = control_input(obj, mu, s)
    % Compute the auxiliary variables
    w = 1+s(2,:).*cos(s(6,:))+s(3,:).*sin(s(6,:));
    s = 1+s(4,:).^2+s(5,:).^2;

    % Compute the control input matrix
    for i = 1:size(s,2)
        B(1,2) = 2*s(1,i)./w(i).*sqrt(s(1,i)./mu);
        B(2,1) = sin(s(6,i)).*sqrt(s(1,i)./mu);
        B(2,2) = sqrt(s(1,i)/mu)./w(i).*((w(i)+1).*cos(s(6,i))+s(2,i));
        B(2,3) = -sqrt(s(1,i)/mu).*s(3,i)./w(i).*(s(4,i).*sin(s(6,i))-s(5,i).*cos(s(6,i)));
        B(3,1) = -cos(s(6,i)).*sqrt(s(1,i)/mu);
        B(3,2) = sqrt(s(1,i)/mu)./w(i).*((w(i)+1).*sin(s(6,i))+s(3,i));
        B(3,3) = sqrt(s(1,i)/mu).*s(2,i)./w(i).*(s(4,i).*sin(s(6,i))-s(5,i).*cos(s(6,i)));
        B(4,3) = sqrt(s(1,i)/mu).*s(i).*cos(s(6,i))./(2*w(i));
        B(5,3) = sqrt(s(1,i)/mu).*s(i).*sin(s(6,i))./(2*w(i));
        B(6,3) = sqrt(s(1,i)/mu)./w(i).*(s(4,i).*sin(s(6,i))-s(5,i).*cos(s(6,i)));
    end
end