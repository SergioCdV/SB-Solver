%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Control function %% 
% Function implementation of the control function as a dynamics residual

function [u] = ControlFunction(obj, params, beta, t0, tf, t, s)
    % Compute the auxiliary variables
    S = 1+s(1,:).*cos(t(1,:))+s(2,:).*sin(t(1,:));

    % Preallocation
    u = zeros(3,size(t,2));                       % Control vector
    u(2,:) = s(9,:) .* (-s(3,:).^3 .* S.^3);      % Y acceleration

    for i = 1:size(t,2)
        u(1,i) = (s(7,i)-(s(1,i)+(1+S(i))*cos(t(1,i))) * u(2,i))^2 + (s(8,i)-(s(2,i)+(1+S(i))*sin(t(1,i))) * u(2,i))^2;

        div = S(i) * sin(t(1,i));
        if (div ~= 0)
            u(1,i) = sqrt(u(1,i) / S(i)^2) * sign( (s(7,i)-(s(1,i)+(1+S(i))*cos(t(1,i)) )*u(2,i)) / div );
        end
    end

    % Compute the complete osculating LVLH quaternion and angular velocity
    omega = zeros(3,size(t,2));

    for i = 1:size(t,2)
        norm_squared = dot(s(4:6,i),s(4:6,i));
        if (norm_squared * norm_squared >= 1)
            s(4:6,i) = -s(4:6,i) / norm_squared;
            norm_squared = dot(s(4:6,i),s(4:6,i));
        end
        B = QuaternionAlgebra.Quat2Matrix([s(4:6,i); -1]);
        omega(:,i) = 4 * B.' * s(10:12,i) / (1 + norm_squared)^2;
    end

    % Normal acceleration
    u(3,:) = sqrt(dot(omega, omega, 1)) .* sign( s(10,:) ./ cos(t(1,:)) );

    % Dimensioning 
    u([1 3],:) = u([1 3],:) .* (s(3,:).^4 .* S.^3);
end