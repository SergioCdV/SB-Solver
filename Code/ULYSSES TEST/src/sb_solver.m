%% Project: 
% Date: 19/05/22

%% Shapse-based optimization %%
% Function to compute the low-thrust orbital transfer using a polynomial
% shape-based approach

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - vector initial_coe, the initial orbital elements 
%         - vector final_coe, the final orbital elements 
%         - scalar K, an initial desired revolutions value 
%         - scalar T, the maximum allowed acceleration 
%         - structure setup, containing the setup of the algorithm in general

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final dV cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - scalar tfapp, the initial estimated time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

    % Setup of the algorithm
    sampling_distribution = 'Chebyshev';     % Sampling grid to be used
    m = 60;                                 % Number of nodes in the grid
    cost = 'Minimum fuel';                  % Cost function to be minimized

    tf = 1; 
    L = 1/6;

    initial = [0 1];
    final = [0 -1];
 
    % Initial guess for the boundary control points
    [tau, W, J, D] = quadrature(3, m, sampling_distribution);
    P0 = repmat(initial.', 1, length(tau));
    
    % Initial guess reshaping
    x0 = reshape(P0, [size(P0,1)*size(P0,2) 1]);
    u0 = zeros(1*length(tau),1);
    x0 = [x0; u0];
    
    % Upper and lower bounds 
    P_lb = [-Inf*ones(size(reshape(P0,[],1))); -Inf*ones(length(tau),1)];
    P_ub = [Inf*ones(size(reshape(P0,[],1))); Inf*ones(length(tau),1)];
    
    % Objective function
    objective = @(x)cost_function(cost, initial.', final.', D, tf, m, tau, W, x);

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

    options = optimoptions('fmincon','Display','off',...
        'Algorithm','sqp','MaxFunEvals',Inf,...
        'OptimalityTolerance', 100*eps, 'StepTolerance', 100*eps);
    
    % Optimisation
    [sol, dV, exitflag, output] = fmincon(objective, x0, A, b, Aeq, beq, P_lb, P_ub, nonlcon, options);
    
    % Solution 
    C = reshape(sol(1:2*(m+1)), [size(P0,1) size(P0,2)]);   % Optimal control points
    u = reshape(sol(end-m:end), [1 size(P0,2)]);    % Optimal control points

    fprintf('Exit flag: %i\n', exitflag)
    if (exitflag ~= 1)
        fprintf("Exit messsage: %s", output.message);
    end

    fprintf("Number of iterations: %i\n", output.iterations);
    fprintf("Number of function evaluations: %i\n", output.funcCount);
    fprintf("Constraint violation: %f \n", output.constrviolation);

    % Compute the optimal analytical solution 
    t = (tau+1)/2;
    uopt = [-(2/3)/L*(1-t(t <= 3*L)/(3*L)) 0*ones(1,length(t(t > 3*L & t <= 1-3*L))) -(2/3)/L*(1-(1-t(t > 1-3*L & t <= 1))/(3*L))];
    
    % Results 
    figure
    plot(tau, C); 
    legend('$x$', '$\dot{x}$')
    xlabel("$t$")
    ylabel("$\mathbf{s}$")
    grid on;

    figure
    hold on
    plot(tau, uopt);
    scatter(tau, u)
    xlabel("$t$")
    ylabel("$\mathbf{u}$")
    grid on;

