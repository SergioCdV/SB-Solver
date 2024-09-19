%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)  

    Tmax = 1.00 * params(2);
    Tmin = 0.97 * params(2);

    % Inequality constraints 
    c = [   
            
            -s(3,:) ...                                             % Positive definiteness of the angular momentum
            +dot(s(4:7,:), s(4:7,:), 1) - 1.1 ...                  % Quaternion norm constraint
            -dot(s(4:7,:), s(4:7,:), 1) + 0.9 ...                  % Quaternion norm constraint
            +dot(u(1:3,:), u(1:3,:), 1) - (Tmax.^2) ...             % Thrust modulation
%             -dot(u(1:2,:), u(1:2,:), 1) + (Tmin.^2) ...             % Thrust modulation
        ].';

    % Equalities 
    beta = atan2(s(2,end), s(1,end));
    theta_f = t(1,end) - beta;

    ceq = [
                dot(s(1:2,end), s(1:2,end), 1) - params(3)^2 ...       % Final eccentricity
                cos(theta_f)-cos( params(5) ) ...                      % Sundman transformation
                sin(theta_f)-sin( params(5) )                          % Sundman transformation
           ];
end