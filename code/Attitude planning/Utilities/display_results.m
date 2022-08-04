%% Project: 
% Date: 31/01/22

%% Display results %%
% Function to display the results of the optimization

% Inputs: - scalar exitflag, the optimisation exitflag
%         - structure output, containing information about the optimisation process 
%         - string cost, the function to be minimized
%         - scalar tfapp, the initial estimated time of flight
%         - scalar tf, the final computed time of flight 
%         - scalar dV, the final optimal cost

function display_results(exitflag, output, cost, tfapp, tf, dV)
    % Constants
    sec2hours = 1/3600;

    % Print the results of the optimisation
    fprintf('Exit flag: %i\n', exitflag)
    if (exitflag ~= 1)
        fprintf("Exit messsage: %s", output.message);
    end

    fprintf("Number of iterations: %i\n", output.iterations);
    fprintf("Number of function evaluations: %i\n", output.funcCount);
    fprintf("Constraint violation: %f \n", output.constrviolation);

    % Time of flight results
    fprintf("Initial estimation of flight time: %0.2f hours\n", tfapp*sec2hours);
    fprintf("Final calculation of flight time: %0.2f hours\n", tf*sec2hours);

    % Cost results
    switch (cost)
        case 'Minimum energy'
            fprintf("Final cost: %0.2f rad/s\n\n", dV);
        case 'Minimum power'
            fprintf("Final cost: %0.2f J\n\n", dV);
        case 'Minimum time'
            fprintf("Final cost: %0.2f hours\n\n", tf*sec2hours);
    end
end