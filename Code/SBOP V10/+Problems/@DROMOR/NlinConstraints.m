%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)    
    % Inequality constraints
    ct = dot(u,u,1)-params(2).^2;           % Thrust modulation
    hpos = -s(3,:);                         % Positive definiteness of the angular momentum
    MRP = dot(s(4:6,:),s(4:6,:),1)-1;       % MRPs norm constraint
    
    c = [ 
            hpos ... 
            MRP ...
            ct
        ].';

    % Equalities 
    beta = atan2(s(2,end), s(1,end));
    theta_f = t(1,end) - beta;

    ceq = [
            dot(s(1:2,end),s(1:2,end),1)-params(3)^2 ...         % Final eccentricity
%             cos(theta_f)-cos(params(5)) ...                      % Sundman transformation
%             sin(theta_f)-sin(params(5))                          % Sundman transformation
           ].';
end