%% Osculating energy %%
% Function to compute the osculating Hamiltonian in the extended phase
% space, for different regularization functions

% Inputs: - scalar mu, the gravitational constant of the problem
%         - array s, the state vector in KS space
%         - string sundman_n, indicating the regularization to be used

% Outputs: - array E, the associated osculating energy

function [E, alpha] = OscEnergy(mu, s, sundman_n)
    % Sanity check and default values 
    if ~exist("sundman_n", "var")
        sundman_n = "1"; 
    end

    % Compute the energy 
    switch sundman_n
        case "1"
            r = dot(s(1:4,:), s(1:4,:), 1);                                 % Radius vector
            E = (2 * dot(s(5:8,:), s(5:8,:), 1) - mu) ./ r;                 % Energy value
            alpha = - E * (2 / mu);                                         % Inverse of the energy

        case "Ecc"
            r = dot(s(1:4,:), s(1:4,:), 1);                                 % Radius vector
            alpha = 2 ./ r ./ (1 + 4 ./ r .* dot(s(5:8,:), s(5:8,:), 1));   % Inverse of the negative energy
            E = - mu / 2 * alpha;                                           % Energy value   

        case "Nu"
            r = dot(s(1:4,:), s(1:4,:), 1);                                 % Radius vector
            E = (2 ./ r.^4 .* dot(s(5:8,:), s(5:8,:), 1) - mu ./ r);        % Energy value
            alpha = - E * (2 / mu);                                         % Inverse of the energy

        case "Cart"
            r = sqrt(dot(s(1:3,:), s(1:3,:), 1));                           % Radius vector
            E = (0.5 * dot(s(4:6,:), s(4:6,:), 1) - mu ./ r);               % Energy value
            alpha = - E * (2 / mu);                                         % Inverse of the energy

        otherwise
            error("The selected regularization does not have an osculatign energy equation... Aborting");
    end

end