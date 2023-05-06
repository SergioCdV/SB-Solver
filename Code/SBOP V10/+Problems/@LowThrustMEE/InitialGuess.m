%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Initial guess function %% 
% Function implementation of the a warming up initial guess if available

function [beta, t0, tf] = InitialGuess(obj, params, initial, final)
    % Constants
    t0 = initial(2);

    % Approximation of the number of revolutions
    dtheta = final(2)-initial(2);
    if (dtheta < 0)
        dtheta = dtheta + 2*pi; 
    end
    
    Napp = ceil( (dtheta+tf*0.5*(initial(4)+final(4)) ) / (2*pi) );
    if (Napp <= 0)
        Napp = 1;
    end 
    
    % New initial TOF
    tf = tf * Napp;
    beta = final(2)+2*pi*Napp;
end