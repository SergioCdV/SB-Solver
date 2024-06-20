%% Project: 
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer using a polynomial
% shape-based approach

% Inputs: - class Problem, defining the problem of interest

% Outputs: - array C, the final state evolution matrix
%          - scalar cost, the final cost of the optimization
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar t0, the initial time of flight 
%          - scalar tf, the final time of flight 
%          - vector t, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [C, cost, u, t0, tf, t, exitflag, output, P] = solve(obj, Problem)
    % Last checks 
    if ((obj.NumNodes + 1) * (Problem.StateDim + Problem.ControlDim) < Problem.StateDim * (max(obj.PolOrder) + 1))
        warning('The problem is over-constrained. Refining the grid...');
        obj.NumNodes = Problem.StateDim * (max(obj.PolOrder) + 1) / (Problem.StateDim + Problem.ControlDim) - 1; 
    end

    if (length(obj.PolOrder) ~= Problem.StateDim)
        warning('The input polynomial order vector mismatches the state dimension...'); 
        obj.PolOrder = [obj.PolOrder min(obj.PolOrder)*ones(1,Problem.StateDim-length(obj.PolOrder))].';
    elseif (size(obj.PolOrder,1) ~= Problem.StateDim)
        obj.PolOrder = obj.PolOrder.';
    end 

    % Setup of the algorithm
    n = obj.PolOrder;                     % Order in the approximation of the state vector
    basis = obj.Basis;                    % Polynomial basis to be used 
    L = Problem.DerDeg;                   % Highest derivative in the dynamics
 
    % Initial guess for the boundary control points
    mapp = 300;   
    Grid = obj.gridding(mapp);

    obj.PolOrder = 3 * ones(size(n));
    B = obj.state_basis(L, obj.PolOrder, basis, Grid.tau);
    [betaapp, t0app, tfapp, ~, Capp] = obj.initial_approximation(Problem, B, Grid); 
    obj.PolOrder = n;
    
    if (obj.InitialGuessFlag)
        P0 = obj.P0;
    else
        % Initial fitting for n+1 control points
        [P0, ~] = obj.initial_fitting(Problem, Grid, Capp);
    end
    
    % Quadrature definition
    Grid = obj.gridding();

    % Final state basis
    [B, CB] = obj.state_basis(L, n, basis, Grid.tau);

    % Initial guess reshaping
    x0 = reshape(P0, size(P0,1) * size(P0,2), []);
    x0 = [x0; t0app; tfapp; betaapp];
        
    % Objective function
%     objective = @(x)obj.cost_function(Problem, B, Grid, x);

    % Non-linear constraints
%     nonlcon = @(x)obj.constraints(Problem, B, B, Grid, x);

    % Declare the global variables 
    global x_last obj_all c_all ceq_all; 
    x_last = x0;
    [obj_all, c_all, ceq_all] = obj.ComputeAll(Problem, B, CB, Grid, x_last);
    ComputeAll = @(x)obj.ComputeAll(Problem, B, CB, Grid, x);

    % Upper and lower bounds 
    [P_lb, P_ub] = obj.opt_bounds(Problem, n, size(betaapp,1));

    % Linear constraints
    [A, b, Aeq, beq] = Problem.LinConstraints(Problem.Params, betaapp, P0);

    % Modification of fmincon optimisation options and parameters (according to the details in the paper)
    options = optimoptions('fmincon', 'TolCon', 1e-6, 'Display', 'off', 'Algorithm', 'sqp', 'ScaleProblem', true);
    options.MaxFunctionEvaluations = obj.maxFunctionEvaluations;
    options.MaxIterations = obj.maxIter;
    
    % Optimisation
    [sol, cost, exitflag, output] = fmincon(@(x)obj_fun(x, ComputeAll), x0, A, b, Aeq, beq, P_lb, P_ub, @(x)non_fun(x, ComputeAll), options);
    
    % Solution 
    StateCard = (max(n)+1) * Problem.StateDim;                               % Cardinal of the state modes
    P = reshape(sol(1:StateCard), Problem.StateDim, []);                     % Control points
    t0 = sol(StateCard+1);                                                   % Initial independent variable value
    tf = sol(StateCard+2);                                                   % Final independent variable value
    beta = sol(StateCard+3:end);                                             % Extra optimization parameters

    [t(1,:), t(2,:)] = Grid.Domain(t0, tf, Grid.tau);                        % Original time independent variable
    
    % Final control points imposing boundary conditions
    P = obj.boundary_conditions(Problem, beta, t0, tf, t, B, P);

    % Final state evolution
    C = obj.evaluate_state(n, L, P, B);

    % Normalization with respect to the independent variable
    m = Problem.StateDim;
    for i = 1:L
        C(1+m*i:m*(i+1),:) = C(1+m*i:m*(i+1),:) ./ ( t(2,:).^i );     
    end

    u = Problem.ControlFunction(Problem.Params, beta, t0, tf, t, C);    % Control function
    t = t(1,:);
    
    % Results 
    % obj.display_results(exitflag, cost, output);
end

%% Auxiliary functions 
% Objective function
function [obj] = obj_fun(x, ComputeAll)
    % Global function 
    global x_last obj_all c_all ceq_all;

    % Check if the function has changed 
    if any(x ~= x_last)
        [obj_all, c_all, ceq_all] = ComputeAll(x);
        x_last = x;
    end

    % Set the new objective function 
    obj = obj_all;
end

% Constraint function
function [c, ceq] = non_fun(x, ComputeAll)
    % Global function 
    global x_last obj_all c_all ceq_all;

    % Check if the function has changed 
    if any(x ~= x_last)
        [obj_all, c_all, ceq_all] = ComputeAll(x);
    end

    % Set the new constraint function 
    c = c_all;
    ceq = ceq_all;
end