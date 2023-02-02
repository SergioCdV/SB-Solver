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
    w = 1+C(2,:).*cos(C(6,:))+C(3,:).*sin(C(6,:));
    s = 1+C(4,:).^2+C(5,:).^2;

    % Linear terms of the equations of motion
    f = tf*[zeros(5,size(C,2)); sqrt(mu.*C(1,:)).*(w./C(1,:)).^2];                  % Acceleration vector
    a = C(7:12,:);                                                                  % Inertial acceleration field
    dv = C(7:12,:);                                                                 % Inertial velocity field

    % Compute the control vector as a dynamics residual
    u = zeros(3,size(C,2));

    % Tangential component
    alpha = 2*C(1,:)./w.*sqrt(C(1,:)/mu); 
    u(2,:) = (a(1,:)-f(1,:))./alpha;

    % Normal component
    beta = sqrt(C(1,:)/mu).*s./(2*w);
    u(3,:) = sqrt(a(4,:).^2+a(5,:).^2)./beta;
    u(3,:) = u(3,:).*sign(a(4,:))./sign(cos(C(6,:)));

    % Radial component
    delta = sqrt(C(1,:)/mu);
    for i = 1:size(C,2)
        B = control_input(mu, C(:,i));
        u(1,i) = (a(2,i)-B(2,2:3)*u(2:3,i))^2+(a(3,i)-B(3,2:3)*u(2:3,i))^2;
    end
    u(1,:) = sqrt(u(1,:))./delta;
end


