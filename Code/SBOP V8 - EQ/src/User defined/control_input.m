%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 22/11/2022

%% Control input matrix %% 
% Function to transform classical orbital elements into Cartesian or
% viceversa

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the orbital equinoctial state vector to be transformed 

% Outputs: - matrix B, the control input matrix

function [B] = control_input(mu, C)
    % Compute the auxiliary variables
    w = 1+C(2,:).*cos(C(6,:))+C(3,:).*cos(C(6,:));
    s = 1+C(4,:).^2+C(5,:).^2;

    % Compute the control input matrix
    for i = 1:size(C,2)
        B(1,2) = 2*C(1,i)./w(i).*sqrt(C(1,i)./mu);
        B(2,1) = sin(C(6,i)).*sqrt(C(1,i)./mu);
        B(2,2) = sqrt(C(1,i)/mu)./w(i).*((w(i)+1).*cos(C(6,i))+C(2,i));
        B(2,3) = -sqrt(C(1,i)/mu).*C(3,i)./w(i).*(C(4,i).*sin(C(6,i))-C(5,i).*cos(C(6,i)));
        B(3,1) = -cos(C(6,i)).*sqrt(C(1,i)/mu);
        B(3,2) = sqrt(C(1,i)/mu)./w(i).*((w(i)+1).*sin(C(6,i))+C(3,i));
        B(3,3) = sqrt(C(1,i)/mu).*C(2,i)./w(i).*(C(4,i).*sin(C(6,i))-C(5,i).*cos(C(6,i)));
        B(4,3) = sqrt(C(1,i)/mu).*s(i).*cos(C(6,i))./(2*w(i));
        B(5,3) = sqrt(C(1,i)/mu).*s(i).*sin(C(6,i))./(2*w(i));
        B(6,3) = sqrt(C(1,i)/mu)./w(i).*(C(4,i).*sin(C(6,i))-C(5,i).*cos(C(6,i)));
    end
end