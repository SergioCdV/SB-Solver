%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 30/01/2022

%% Compute All %% 
% Function to compute the residual vector of the constraints of the
% problem, together with the objective function

% Inputs: - class Problem, defining the problem at hands
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - vector tau, the normalized independent variable
%         - vector x, the degree of freedom to be optimized

% Outputs: - objective function r
%          - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [r, c, ceq] = ComputeAll(obj, Problem, B, CB, Grid, x)
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
    s = obj.evaluate_state(n, L, P, CB);                                % State evolution
    
    % Normalization
    for i = 1:L
        s(1+m*i:m*(i+1),:) = s(1+m*i:m*(i+1),:) ./ ( t(2,:).^i ); 
    end

    u = Problem.ControlFunction(Problem.Params, beta, t0, tf, t, s);    % Control function

    % Equalities 
    [c, ceq] = Problem.NlinConstraints(Problem.Params, beta, t0, tf, [t(1,:); t(2,:) .* Grid.W], s, u);

%     % Relaxation 
%     N = length(Grid.tau);
%     delta = 0;%(N-1)^(1/2-L);
%     c = [c-delta; -c-delta; ceq-delta; -ceq-delta];
%     ceq = [];

    % Evaluate the cost function (Lagrange and Mayer terms)
    [M, L] = Problem.CostFunction(Problem.Params, beta, t0, tf, t(1,:), s, u); 

    if (isempty(Grid.W))
        r = M + trapz(t(2,:).*Grid.tau, L);
    else
        r = M + dot(t(2,:).*Grid.W, L);
    end
end