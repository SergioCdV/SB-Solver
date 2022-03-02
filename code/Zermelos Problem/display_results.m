%% Project: 
% Date: 31/01/22

%% Display results %%
% Function to display the results of the optimization

% Inputs: - array P0, the initial control points in the trajectory Bernstein
%           approximation
%         - array P, the final control points in the trajectory Bernstein
%           approximation
%         - array B, the Bernstein polynomial basis used in the
%           approximation 
%         - scalar exitflag, the optimisation exitflag
%         - structure output, containing information about the optimisation process 
%         - scalar tfapp, the initial estimated time of flight 
%         - scalar r0, the dimensionalising characteristic distance 
%         - scalar n, the order of the approximation 

function display_results(P0, P, B, m, exitflag, output, tfapp, r0, n)
    % Constants
    days2sec = 86400;

    % Print the results of the optimisation
    fprintf('Exit flag: %i\n', exitflag)
    if (exitflag ~= 1)
        fprintf("Exit messsage: %s", output.message);
    end

    fprintf("Number of iterations: %i\n", output.iterations);
    fprintf("Number of function evaluations: %i\n", output.funcCount);
    fprintf("Constraint violation: %f \n", output.constrviolation);

    % Time of flight results
    tf(1) = flight_time(P0, B, m, tfapp, r0, n);
    tf(2) = flight_time(P, B, m, tfapp, r0, n);
    fprintf("Initial estimation of flight time: %0.2f days\n", tfapp/days2sec);
    fprintf("Initial calculation of flight time: %0.2f days\n", tf(1)/days2sec);
    fprintf("Final calculation of flight time: %0.2f days\n\n", tf(2)/days2sec);
end