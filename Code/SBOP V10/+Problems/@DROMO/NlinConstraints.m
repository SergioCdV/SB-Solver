%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Constraints function %% 
% Function implementation of the path and boundary constraints functions

function [c, ceq] = NlinConstraints(obj, params, beta, t0, tf, t, s, u)    
    % Inequality constraints
    ct = dot(u,u,1)-params(2).^2;   % Thrust modulation

    % Equalities (Sundman transformation)
    rvd_cns = OrbitalDynamics.dromo2coe([s(1:7,end); t(end)]);
    ceq = [dot(s(4:7,:),s(4:7,:),1)-1 ...
           rvd_cns(2)-params(4) ...
           cos(rvd_cns(end))-cos(params(5)) ...
           sin(rvd_cns(end))-sin(params(5))
           ].';

    c = [-s(3,:) ct].';
end