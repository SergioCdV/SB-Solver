%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Compute the auxiliary variables
    S = 1+s(1,:).*cos(t(1,:))+s(2,:).*sin(t(1,:));

    % Preallocation
    u = zeros(3,size(t,2));                        % Control vector
    u(2,:) = s(10,:) .* (-s(3,:).^3 .* S.^3);      % Y acceleration

    for i = 1:size(t,2)
        u(1,i) = (s(8,i)-(s(1,i)+(1+S(i))*cos(t(1,i))) * u(2,i))^2 + (s(9,i)-(s(2,i)+(1+S(i))*sin(t(1,i))) * u(2,i))^2;
        u(1,i) = sqrt(u(1,i) / S(i)^2) * sign( (s(8,i)-(s(1,i)+(1+S(i))*cos(t(1,i)) )*u(2,i)) / (S(i) * sin(t(1,i))) );
    end

    % Compute the complete osculating LVLH quaternion and angular velocity
    epsilon = s(4:7,:);

    omega = 2 * [epsilon(4,:).*s(11,:)+epsilon(3,:).*s(12,:)-epsilon(2,:).*s(13,:)-epsilon(1,:).*s(14,:); ...
                -epsilon(3,:).*s(11,:)+epsilon(4,:).*s(12,:)+epsilon(1,:).*s(13,:)-epsilon(2,:).*s(14,:); ...
                 epsilon(2,:).*s(11,:)-epsilon(1,:).*s(12,:)+epsilon(4,:).*s(13,:)-epsilon(3,:).*s(14,:)];

    % Normal acceleration
%     u(3,:) = sqrt(dot(omega, omega, 1)) .* sign( s(11,:) ./ cos(t(1,:)) );
    u(3,:) = omega(1,:) ./ cos(t(1,:));

    % Dimensioning 
    u([1 3],:) = u([1 3],:) .* (s(3,:).^4 .* S.^3);
end