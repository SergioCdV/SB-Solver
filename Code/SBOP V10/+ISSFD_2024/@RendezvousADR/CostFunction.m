%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Differential time law 
    rho = 1 + params(5) * cos(t(1,:));      % Transformation parameter
    k = params(2)^2/params(4)^3;
    Omega = k .* rho.^2;                    % True anomaly angular velocity [rad/s]

    % Mayer and Lagrange terms
    if ( length(params) <= 35 + 4 * (params(35)+1) )
        q = params(35 + 4 * params(35) : 35 + 4 * ( params(35) + 1) ).';
        sigma = QuaternionAlgebra.MPR2Quat(1, 1, q, false);
        dsigma = 0.25 * QuaternionAlgebra.Quat2Matrix([sigma; -1]) * params(32:34).';

        beta_ref = [
                    zeros(3,1); ...
                    sigma; ...
                    zeros(3,1); ...
                    dsigma
                ]; 
        diff = zeros(3,1);
    else
        beta_ref = [
                    params(35 + 4 * (params(35)+1) + 1 : 35 + 4 * (params(35)+1) + 9).';
                    ];
        beta_ref = [beta_ref(1:6); zeros(3,1); beta(7:9)];
        diff = s(1:12,end) - beta_ref;
    end
    
    M = dot(diff, diff);

    L = dot(u(1:3,:), u(1:3,:), 1) + dot(u(4:6,:), u(4:6,:), 1);     
    L = L .* Omega;
end