%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - class Problem, defining the problem at hands
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - vector tau, the normalized independent variable
%         - vector x, the degree of freedom to be optimized

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(obj, Problem, B, Grid, x)
    % Optimization variables
    L = Problem.DerDeg;                                                 % Maximum derivative degree
    m = Problem.StateDim;                                               % State dimension
    n = obj.PolOrder;                                                   % Order of the polynomial approximation
    StateCard = (max(n)+1) * m;                                         % Cardinal of the state modes
    P = reshape(x(1:StateCard), m, []);                                 % Control points
    t0 = x(StateCard+1);                                                % Initial independent variable value
    tf = x(StateCard+2);                                                % Final independent variable value
    beta = x(StateCard+3:end);                                          % Extra optimization parameters

    % Evaluate the boundary conditions and the state evolution
    [t(1,:), t(2,:)] = Grid.Domain(t0, tf, Grid.tau);                   % Original time independent variable
    P = obj.boundary_conditions(Problem, beta, t0, tf, t, B, P);        % Boundary conditions control points
    s = obj.evaluate_state(n, L, P, B);                                 % State evolution
    
    % Normalization
    for i = 1:L
        s(1+m*i:m*(i+1),:) = s(1+m*i:m*(i+1),:) ./ ( t(2,:).^i );     
    end

    u = Problem.ControlFunction(Problem.Params, beta, t0, tf, t, s);    % Control function

    % Equalities 
    [c, ceq] = Problem.NlinConstraints(Problem.Params, beta, t0, tf, t(1,:), s, u);

%     % Relaxation 
%     N = length(Grid.tau);
%     delta = 0;%(N-1)^(1/2-L);
%     c = [c-delta; -c-delta; ceq-delta; -ceq-delta];
%     ceq = [];
end