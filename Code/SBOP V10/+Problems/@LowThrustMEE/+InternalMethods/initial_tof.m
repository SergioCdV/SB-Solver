%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 16/04/22

%% Initial TOF %%
% Function to estimate the initial time of flight

% Inputs: - scalar mu, the gravitational parameter of the system  
%         - scalar T, the maximum acceleration allowed for the
%           spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the final boundary conditions of the
%           trajectory

% Outputs:- scalar tfapp, the initial approximation of the time of flight

function [tfapp] = initial_tof(obj, mu, T, initial, final)
    % Approximation of the time of flight via the Vis-Viva theorem
    r0 = norm(initial([1 3]));              % Initial radius
    v0 = norm(initial(4:6));                % Initial velocity
    rf = norm(final([1 3]));                % Final radius
    vf = norm(final(4:6));                  % Final velocity
    a = (r0+rf)/2;                          % Transfer semimajor axis
    dE = (vf^2-v0^2)/2-mu*(1/rf-1/r0);      % Change of energy
    tfapp = 2*abs(dE)/(T*sqrt(mu/a));       % Time of flight
end