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

function [c, ceq] = constraints(Problem, B, basis, domain_mapping, tau, x)
    % Optimization variables
    StateCard = (max(Problem.PolOrder)+1) * Problem.StateDim;                % Cardinal of the state modes
    P = reshape(x(1:StateCard), Problem.StateDim, []);                       % Control points
    t0 = x(StateCard+1);                                                     % Initial independent variable value
    tf = x(StateCard+2);                                                     % Final independent variable value
    beta = x(StateCard+3:end);                                               % Extra optimization parameters

    % Evaluate the boundary conditions and the state evolution
    P = boundary_conditions(Problem, beta, t0, tf, B, basis, P);             % Boundary conditions control points
    s = evaluate_state(P, B, Problem.PolOrder);                              % State evolution

    % Evaluate the control function 
    [t, J] = feval(domain_mapping, t0, tf, tau);                             % Original time independent variable and Jacobian of the transformation
    u = Problem.ControlFunction(Problem.ControlDim, beta, t0, tf, t, s);     % Control function
    u = J * u;

    % Equalities 
    [c, ceq] = Problem.Constraints(beta, t0, tf, tau, s, u);
end