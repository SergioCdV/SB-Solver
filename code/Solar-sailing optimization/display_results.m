function [] = display_results(P0, P, B, m, exitflag, output, tfapp, r0, n)

% DISPLAY_RESULTS prints to the terminal a short summary of the
% optimisation results.

seconds_per_day = 60*60*24;

fprintf('Exit flag: %i\n', exitflag)
if (exitflag ~= 1)
    fprintf("Exit messsage: %s\n", output.message)
end



fprintf("Number of iterations: %i\n", output.iterations)
fprintf("Number of function evaluations: %i\n", output.funcCount)
fprintf("Constraint violation: %f \n", output.constrviolation)


% see variation of flight time calculations

fprintf("Initial estimation of flight time: %0.2f days\n", tfapp/seconds_per_day)
fprintf("Initial calculation of flight time: %0.2f days\n", flight_time(P0, B, m, tfapp, r0, n)/seconds_per_day)
fprintf("Final calculation of flight time: %0.2f days\n\n", flight_time(P, B, m, tfapp, r0, n)/seconds_per_day)





end