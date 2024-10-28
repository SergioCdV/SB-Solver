%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)  

    Tmax = 1.00 * params(2);
%     Tmin = 0.97 * params(2);

    zeta = atan2( beta(2), beta(1) );
    alpha = t(1,end) - zeta;
    theta_f = atan2( sin(alpha), cos(alpha) );

    % Inequality constraints 
    c = [   
%             -s(3,:)...
            +dot(u(1:3,:), u(1:3,:), 1) - (Tmax.^2) ...             % Thrust modulation
            +dot(s(4:7,:), s(4:7,:), 1) - 1 ...                     % Quaternion norm constraint
            -dot(s(4:7,:), s(4:7,:), 1) + 0.98 ...                  % Quaternion norm constraint
%             -dot(u(1:3,:), u(1:3,:), 1) + (Tmin.^2) ...           % Thrust modulation
        ].';

    % Equalities 
%     l = 1 + s(1,:) .* cos(t(1,:)) + s(2,:) .* sin(t(1,:));
%     gamma = s(3,:).^3 .* l.^2;
    ceq = [
%                 50 - trapz(t(1,:), 1./gamma)
                +dot(beta, beta) - params(3)^2 ...                      % Final eccentricity
                +cos( theta_f )-cos( params(5) ) ...                    % Final osculating true anomaly
                +sin( theta_f )-sin( params(5) )  ...                   % Final osculating true anomaly
           ];
end