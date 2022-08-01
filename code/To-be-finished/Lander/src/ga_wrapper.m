%% Project: 
% Date: 19/05/22

%% GA wrapper %%
% Function to compute the optimal basis and grid using metaheuristic GA and
% the associated optimal results

% Inputs: - structure system, containing the physical information of the
%           2BP of interest
%         - vector initial_coe, the initial orbital elements 
%         - vector final_coe, the final orbital elements 
%         - scalar K, an initial desired revolutions value 
%         - scalar T, the maximum allowed acceleration
%         - structure setup, containing the setup of the figures

% Outputs: - array C, the final state evolution matrix
%          - scalar dV, the final dV cost of the transfer 
%          - array u, a 3xm matrix with the control input evolution  
%          - scalar tf, the final time of flight 
%          - scalar tfapp, the initial estimated time of flight 
%          - vector tau, the time sampling points final distribution
%          - exitflag, the output state of the optimization process 
%          - structure output, containing information on the final state of
%            the optimization process

function [sol, C, dV, u, tf, tfapp, tau, exitflag, output] = ga_wrapper(system, initial_coe, final_coe, K, T, setup)
    % Genetic algorithm setup
    dof = 4;               % Number of DOF of the optimization process
    PopSize = 10;         % Population size for each generation
    MaxGenerations = 10;   % Maximum number of generations for the evolutionary algorithm
            
    options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);

    intcon = 1:dof;        % All DOF are integer

    A = []; 
    b = []; 
    Aeq = [];
    beq = [];
    nonlcon = [];

    % Lower and upper bounds 
    lb = [30 1 1 5];
    ub = [100 7 6 15];

    % Metaheuristic selection 
    sol = ga(@(x)gasp_opti(system, initial_coe, final_coe, K, T, x, setup), dof, A, b, Aeq, beq, lb, ub, nonlcon, intcon, options);

    % Final results 
    m = sol(1);                         % Number of sampling points 
    sampling_distribution = sol(2);     % Final sampling distribution 
    switch (sampling_distribution)
        case 1
            sampling_distribution = 'Linear';
        case 2 
            sampling_distribution = 'Normal';
        case 3 
            sampling_distribution = 'Random';
        case 4 
            sampling_distribution = 'Legendre';
        case 5 
            sampling_distribution = 'Chebyshev';
        case 6
            sampling_distribution = 'Regularized';
    end

    basis = sol(3);                     % Final sampling distribution 
    switch (basis)
        case 1
        case 2 
    end
    n = sol(3);                         % Final approximation order

    sol.Basis = basis; 
    sol.Order = n; 
    sol.Points = sampling_distribution; 
    sol.NumPoints = m; 
    
    % Optimal results
    [C, dV, u, tf, tfapp, tau, exitflag, output] = spaed_optimization(system, initial_coe, final_coe, K, T, m, sampling_distribution, basis, n, setup);
end
 

%% Auxiliary functions 
% Function to determine the best configuration 
function [cost] = gasp_opti(system, initial_coe, final_coe, K, T, x, setup)
    % Setup 
    m = x(1);                         % Number of sampling points 
    sampling_distribution = x(2);     % Final sampling distribution 
    switch (sampling_distribution)
        case 1
            sampling_distribution = 'Linear';
        case 2 
            sampling_distribution = 'Normal';
        case 3 
            sampling_distribution = 'Random';
        case 4 
            sampling_distribution = 'Legendre';
        case 5 
            sampling_distribution = 'Chebyshev';
        case 6
            sampling_distribution = 'Regularized';
        case 7
            sampling_distribution = 'Laguerre';
    end

    basis = x(3);                     % Final sampling distribution 
    switch (basis)
        case 1
            basis = 'Bernstein';
        case 2 
            basis = 'Orthogonal Bernstein';
        case 3
            basis = 'Legendre';
        case 4
            basis = 'Chebyshev';
        case 5
            basis = 'Hermite';
        case 6
            basis = 'Laguerre';
    end
    n = x(4);                         % Final approximation order

    % Optimization
    [~, dV, ~, tf, ~, ~, exitflag, ~] = spaed_optimization(system, initial_coe, final_coe, K, T, m, sampling_distribution, basis, n, setup);

    % Final cost 
    cost = dV*tf/((exitflag == 1) || (exitflag == 2) || (exitflag == -3));
end
