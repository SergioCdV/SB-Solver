%% Project: 
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer using a polynomial
% shape-based approach

% Inputs: - vector initial, the initial state 
%         - vector final, the final state
%         - scalar tf, the time of flight
%         - scalar L, a constrain on the X dimension
%         - structure setup, containing the setup of the algorithm in general

% Outputs: - array C, the final state evolution matrix
%          - scalar cost, the minimum cost function value
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, cost, u, tf, tau, exitflag, output] = sb_solver(initial, final, tf, L, setup)
    % Setup of the algorithm
    sampling_distribution = setup.grid;     % Sampling grid to be used
    m = setup.nodes;                        % Number of nodes in the grid
 
    % Initial guess for the boundary control points
    [tau, W, ~, D] = quadrature(3, m, sampling_distribution);
    P0 = repmat(initial.', 1, length(tau));
    
    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    u0 = zeros(1*length(tau),1);
    x0 = [x0; u0];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(size(reshape(P0,[],1))); -Inf*ones(length(tau),1)];
    P_ub = [Inf*ones(size(reshape(P0,[],1))); Inf*ones(length(tau),1)];
    
    % Objective function
    objective = @(x)cost_function(tf, m, tau, W, x);

    % Non-linear constraints
    nonlcon = @(x)constraints(initial.', final.', L, tf, D, m, tau, x);
    
    % Linear constraints and inequalities
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp');
    options.MaxFunctionEvaluations = 1e6;

    options = optimoptions('fmincon','Display','off', 'Algorithm', 'sqp', 'MaxFunEvals', Inf, 'OptimalityTolerance', 100*eps, 'StepTolerance', 100*eps);
    
    % Optimisation
    [sol, cost, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    C = reshape(sol(1:2*(m+1)), [size(P0,1) size(P0,2)]);   % Optimal state evolution
    u = reshape(sol(end-m:end), [1 size(P0,2)]);            % Optimal control law

    tau = 0.5*(tau+1);

    fprintf('Exit flag: %i\n', exitflag)
    if (exitflag ~= 1)
        fprintf("Exit messsage: %s", output.message);
    end

    fprintf("Number of iterations: %i\n", output.iterations);
    fprintf("Number of function evaluations: %i\n", output.funcCount);
    fprintf("Constraint violation: %f \n", output.constrviolation);
end
