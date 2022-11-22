%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/04/22

%% Acceleration control %%
% Function to compute the acceleration vector norm from the nominal
% trajectory

% Inputs: - scalar mu, the gravitational parameter of the system
%         - array C, the 9xm state vector 
%         - scalar tf, the final time of flight

% Outputs: - vector u, the nondimensional 3xm control vector
%          - vector dv, the inertial velocity field 
%          - vector f, the dynamics vector field

function [u, dv, f] = acceleration_control(mu, C, tf)
    % Compute the auxiliary variables
    w = 1+C(2,:).*cos(C(6,:))+C(3,:).*cos(C(6,:));
    s = 1+C(4,:).^2+C(5,:).^2;

    % Control input matrix 
    B = zeros(6,3); 

    % Linear terms of the equations of motion
    f = tf*[zeros(5,size(C,2)); sqrt(mu.*C(1,:)).*(w./C(1,:)).^2];                  % Acceleration vector
    a = C(7:12,:);                                                                  % Inertial acceleration field
    dv = C(7:12,:);                                                                 % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = zeros(3,size(C,2));

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
        u(:,i) = B\(a(:,i)-f(:,i));
    end
end
