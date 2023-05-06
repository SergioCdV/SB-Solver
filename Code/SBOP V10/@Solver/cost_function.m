%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Cost function %%
% Function to compute the cost function to be minimized

% Inputs: - class Problem, defining the problem at hands
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - function handle domain_mapping, to map the independent variable
%           domain to the original one
%         - vector tau, the vector of collocation points
%         - vector W, the quadrature weights
%         - vector x, the degree of freedom to be optimized
%         - string dynamics, the independent variable parametrization to be
%           used

% Outputs: - scalar r, the cost index to be optimized

function [r] = cost_function(obj, Problem, B, Grid, x)
    % Optimization variables
    L = Problem.DerDeg;                                                 % Maximum derivative degree
    m = Problem.StateDim;                                               % State dimension
    n = obj.PolOrder;                                                   % Order of the polynomial approximation
    StateCard = (max(n)+1) * m;                                         % Cardinal of the state modes
    P = reshape(x(1:StateCard), m, []);                                 % Control points
    t0 = x(StateCard+1);                                                % Initial independent variable value
    tf = x(StateCard+2);                                                % Final independent variable value
    beta = x(StateCard+3:end);                                          % Extra optimization parameters
    
    % Evaluate the boundary conditions
    [t(1,:), t(2,:)] = Grid.Domain(t0, tf, Grid.tau);                   % Original time independent variable
    P = obj.boundary_conditions(Problem, beta, t0, tf, t, B, P);        % Boundary conditions control points
    s = obj.evaluate_state(n, L, P, B);                                 % State evolution

    % Normalization
    for i = 1:L
        s(1+m*i:m*(i+1),:) = s(1+m*i:m*(i+1),:) ./ ( t(2,:).^i );     
    end

    u = Problem.ControlFunction(Problem.Params, beta, t0, tf, t, s);    % Control function
        
    % Evaluate the cost function (Lagrange and Mayer terms)
    [M, L] = Problem.CostFunction(Problem.Params, beta, t0, tf, s, u); 

    if (isempty(Grid.W))
        r = M + trapz(t(2,:).*Grid.tau, L);
    else
        r = M + dot(t(2,:).*Grid.W, L);
    end
end