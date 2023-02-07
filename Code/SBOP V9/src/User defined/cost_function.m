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

function [r] = cost_function(Problem, B, basis, domain_mapping, tau, W, x)
    % Optimization variables
    StateCard = (max(Problem.PolOrder)+1) * Problem.StateDim;           % Cardinal of the state modes
    P = reshape(x(1:StateCard), Problem.StateDim, []);                  % Control points
    t0 = x(StateCard+1);                                                % Initial independent variable value
    tf = x(StateCard+2);                                                % Final independent variable value
    beta = x(StateCard+3:end);                                          % Extra optimization parameters
    
    % Evaluate the boundary conditions
    P = boundary_conditions(Problem, beta, t0, tf, B, basis, P);        % Boundary conditions control points
    s = evaluate_state(P, B, Problem.PolOrder);                         % State evolution

    % Evaluate the control function 
    t = feval(domain_mapping, t0, tf, tau);                                  % Original time independent variable
    u = Problem.ControlFunction(Problem.ControlDim, beta, t0, tf, t, s);     % Control function
        
    % Evaluate the cost function (Lagrange and Mayer terms)
    [L, M] = Problem.evaluate_cost(beta, t0, tf, s, u); 

    if (isempty(W))
        r = M + (tf-t0) * trapz(tau,L);
    else
        r = M + (tf-t0) * dot(W,L);
    end
end
