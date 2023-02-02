%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 31/01/22

%% Initial approximation %%
% Function to estimate the initial time of flight, control points and curve approximation

% Inputs: - scalar mu, the gravitational parameter of the system
%         - vector tau, the collocation points to be used 
%         - vector initial, the initial boundary conditions of the
%           trajectory
%         - vector final, the final boundary conditions of the
%           trajectory
%         - string basis, specifying the polynomial collacation basis

% Outputs: - array Papp, the initial estimation of the boundary control
%            points
%          - initial estimation of the control law
%          - scalar tfapp, the initial initial time of flight

function [Uapp, tfapp] = initial_approximation(mu, tau, tfapp, initial, final, basis)
    % Constants 
    n_init = [3 3 3];

    % Preliminary number of revolutions
    dtheta = final(2)-initial(2);
    if (dtheta < 0)
        dtheta = dtheta + 2*pi; 
    end
    
    Napp = ceil( (dtheta+tfapp*0.5*(initial(4)+final(4)) ) / (2*pi) );
    if (Napp <= 0)
        Napp = 1;
    end 
    
    % New initial TOF
    tfapp = tfapp*Napp;

    % Initial state vector approximation
    Bapp = state_basis(n_init, tau, basis);
    Papp = zeros(length(initial)/2, max(n_init)+1);  
    Papp = boundary_conditions(tfapp, n_init, initial, final, Papp, Bapp, basis);

    N = size(Papp,1);                      % Number of state variables
    Capp = zeros(3*N,size(Bapp{1},2));     % Preallocation for speed

    for i = 1:N
        % State vector fitting
        k = n_init(i)+1;
        for j = 1:3
            Capp(i+N*(j-1),:) = Papp(i,1:n_init(i)+1)*Bapp{i}(1+k*(j-1):k*j,:);
        end
    end

    % Initial estimate of control points (using the non-orthonormal boundary conditions)
    Uapp = Capp(7:9,:)+mu.*Capp(1:3,:)./sqrt(dot(Capp(1:3,:), Capp(1:3,:), 1)).^3;
    Uapp = ones(3,size(Capp,2));
end