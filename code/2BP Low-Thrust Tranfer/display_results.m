%% Project: 
% Date: 04/03/22

%% Display results %%
% Function to display the results of the optimization

% Print the results of the optimisation
fprintf('Exit flag: %i\n', exitflag)
if (exitflag ~= 1)
    fprintf("Exit messsage: %s", output.message);
end

fprintf("Number of iterations: %i\n", output.iterations);
fprintf("Number of function evaluations: %i\n", output.funcCount);
fprintf("Constraint violation: %f \n", output.constrviolation);

% Time of flight results
fprintf("Initial estimation of flight time: %0.2f \n", tfapp);
fprintf("Final calculation of flight time: %0.2f \n", tf);
