%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)    
    % Initial guess
    t0 = params(5);                    % Initial true anomaly [rad]
    tf = params(6);                    % Final true anomaly [rad]

    if ( length(params) <= 35 + 4 * (params(35)+1) )
        q = params(35 + 4 * params(35) : 35 + 4 * ( params(35) + 1) ).';
        sigma = QuaternionAlgebra.MPR2Quat(1, 1, q, false);
        dsigma = 0.25 * QuaternionAlgebra.Quat2Matrix([sigma; -1]) * params(32:34).';

        beta = [
                    zeros(3,1); ...
                    sigma; ...
                    dsigma
                ]; 
    else
        beta = [
                    params(35 + 4 * (params(35)+1) + 1 : 35 + 4 * (params(35)+1) + 9).'; ...
                ];
    end
end