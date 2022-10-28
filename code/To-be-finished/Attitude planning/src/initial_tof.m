%% Project: Shape-based attitude planning %%
% Date: 16/04/22

%% Initial TOF %%
% Function to estimate the initial maneuver time

% Inputs: - matrix I, the gravitational parameter of the system  
%         - scalar T, the maximum acceleration allowed for the
%           spacecraft
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the final boundary conditions of the
%           trajectory

% Outputs:- scalar tfapp, the initial approximation of the time of flight

function [tfapp] = initial_tof(I, T, initial, final)
    % Approximation of the maneuver time via the vis-viva theorem
    q0 = initial(1:4);                             % Initial radius
    qf = final(1:4);                               % Final radius
    qf = quaternion_inverse(qf);
    dq = quaternion_product(qf,q0);                % Quaternion difference

    omega0 = initial(4:6);                         % Initial velocity
    omegaf = final(4:6);                           % Final velocity

    dE = omegaf.'*I*omegaf-omega0.'*I*omega0;      % Change of kinetic energy
    dE = dE + acos(dq(1));                         % Change of potential energy        

    % Maneuver time
    tfapp = abs(dE)/T;
end