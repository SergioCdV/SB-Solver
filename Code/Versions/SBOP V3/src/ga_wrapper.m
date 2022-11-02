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

% Outputs: - structure sol, containing information on the final state of
%            the optimization process
%          - array fval, the functional evaluation of the Pareto front

function [Sol, fval] = ga_wrapper(system, initial_coe, final_coe, K, T, setup)
    % Genetic algorithm setup
    dof = 4;               % Number of DOF of the optimization process
    PopSize = 20;          % Population size for each generation
    MaxGenerations = 20;   % Maximum number of generations for the evolutionary algorithm
            
    options = optimoptions(@ga,'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);

    intcon = 1:dof;        % All DOF are integer

    A = []; 
    b = []; 
    Aeq = [];
    beq = [];
    nonlcon = [];

    % Lower and upper bounds 
    lb = [30 1 1 5];
    ub = [60 4 4 10];

    % Metaheuristic selection 
    [sol, fval] = gamultiobj(@(x)gasp_opti(system, initial_coe, final_coe, K, T, x, setup), dof, A, b, Aeq, beq, lb, ub, nonlcon, intcon, options);

    % Final results 
    for i = 1:size(sol,1)
        m = sol(i,1);                         % Number of sampling points 
        sampling_distribution = sol(i,2);     % Final sampling distribution
        basis = sol(i,3);                     % Final sampling 
        n = sol(i,4);                         % Final approximation order

        % Write to csv 
        file = 'Results\metaheuristic.csv'; 
        if (isempty(dir(file)))
            cHeader = {'m' 'n' 'basis' 'time_distribution'};
            commaHeader = [cHeader; repmat({','}, 1, numel(cHeader))]; 
            commaHeader = commaHeader(:)';
            textHeader = cell2mat(commaHeader);
            
            fid = fopen(file, 'w'); 
            fprintf(fid, '%s\n', textHeader);
            fclose(fid);
        end

        dlmwrite(file, [m n basis sampling_distribution], '-append');
        
        switch (sampling_distribution)
            case 1
                sampling_distribution = 'Linear';
            case 2 
                sampling_distribution = 'Legendre';
            case 3 
                sampling_distribution = 'Chebyshev';
            case 4
                sampling_distribution = 'Regularized';
        end
    
        switch (basis)
            case 1
                basis = 'Bernstein';
            case 2 
                basis = 'Orthogonal Bernstein';
            case 3
                basis = 'Legendre';
            case 4
                basis = 'Chebyshev';
        end
        
        Sol.Basis{i} = basis; 
        Sol.Order(i) = n; 
        Sol.Points{i} = sampling_distribution; 
        Sol.NumPoints(i) = m; 
    end
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
            sampling_distribution = 'Legendre';
        case 3 
            sampling_distribution = 'Chebyshev';
        case 4
            sampling_distribution = 'Regularized';
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
    end
    n = x(4);                         % Final approximation order

    % Optimization
    [~, dV, ~, tf, ~, ~, exitflag, ~] = spaed_optimization(system, initial_coe, final_coe, K, T, m, sampling_distribution, basis, n, setup);

    % Final cost 
    cost = [dV tf/((exitflag == 1) || (exitflag == 2) || (exitflag == -3))];
end
