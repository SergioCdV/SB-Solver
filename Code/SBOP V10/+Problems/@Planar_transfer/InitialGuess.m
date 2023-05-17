%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)
    % Constants
    mu = params(1); 
    T = params(2); 

    % Approximation of the time of flight via the Vis-Viva theorem
%     r0 = norm(initial([1 3]));              % Initial radius
%     v0 = norm(initial(4:6));                % Initial velocity
%     rf = norm(final([1 3]));                % Final radius
%     vf = norm(final(4:6));                  % Final velocity
%     a = (r0+rf)/2;                          % Transfer semimajor axis
%     dE = (vf^2-v0^2)/2-mu*(1/rf-1/r0);      % Change of energy
%     tf = 2*abs(dE)/(T*sqrt(mu/a));          % Time of flight
    t0 = 0;

    % Approximation of the number of revolutions
%     dtheta = final(2)-initial(2);
%     if (dtheta < 0)
%         dtheta = dtheta + 2*pi; 
%     end
%     
%     Napp = ceil( (dtheta+tf*0.5*(initial(4)+final(4)) ) / (2*pi) );
%     if (Napp <= 0)
%         Napp = 1;
%     end 
%     
%     % New initial TOF
    tf = 2;%tf*Napp;
    beta = final(2)+2*pi*2;
end