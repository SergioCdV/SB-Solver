%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - scalar mu, the gravitational parameter of the central body 
%         - scalar T, the maximum acceleration allowed for the spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - string dynamics, the independent variable parametrization to be
%           used
%         - string cost, the cost function to be minimized

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(mu, initial, final, tf, time_free, B, basis, n, tau, x)
    % Extract the optimization variables
    P = reshape(x(1:end-5), [length(n), max(n)+1]);     % Control points
    thetaf = x(end-4);                                  % Final fiber angle
    dE0 = x(end-2);                                     % Initial energy derivative
    dEf = x(end-3);                                     % Final energy derivative
    sf = x(end-1);                                      % Final time of flight 
    T = x(end);                                         % Needed thrust vector

    % Boundary conditions points
    R = [cos(thetaf) 0 0 -sin(thetaf); 0 cos(thetaf) sin(thetaf) 0; 0 -sin(thetaf) cos(thetaf) 0; sin(thetaf) 0 0 cos(thetaf)];

    initial_u = [initial dE0];
    final_u = [final dEf];
    final_u(1:4) = final_u(1:4)*R.';
    final_u(6:9) = final_u(6:9)*R.';

    P = boundary_conditions(sf, n, initial_u, final_u, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P,B,n);

    % Sundman transformation
    r = dot(C(1:4,:),C(1:4,:),1);
    C(5,:) = 2./r./(1+4./r.*dot(C(6:9,:),C(6:9,:),1)/sf^2);

    % Control input 
    [u, ~] = acceleration_control(mu, C, sf);

    % Equalities 
    l = bilinear_function(C(1:4,:),C(6:9,:));
    m = bilinear_function(C(11:14,:),C(1:4,:));
    dE = (mu/8)*sf^2*r.*C(10,:)+dot(C(6:9,:),u,1);

    ceq = [l m];
    
    if (time_free)
        ceq = [ceq tf-sf*trapz(tau, r./sqrt(mu*C(5,:)))];
    end

    % Inequalities
    c = [dot(u,u,1)-sf^2/4*(r.*repmat(T,1,size(u,2))).^2];
end