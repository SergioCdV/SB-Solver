%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 01/02/22

%% Initial fitting %%
% Function to estimate the trajectory approximation

% Inputs: - object Problem, defining the transcription of the problem of interest
%         - cell array B, containing the polynomial support of the state
%           vector
%         - object Grid, the collocation points object to be used 
%         - array s, the initial state vector trajectory

% Outputs: - array P, the estimation of the state points as an array
%          - array C, the estimation of the problem state vector across the
%            grid

function [P, C] = initial_fitting(obj, Problem, Grid, s)
    % Constants 
    n = obj.PolOrder;       % Polynomial order for each state variable 
    L = Problem.DerDeg;     % Order of the dynamics (maximum derivative order)
    basis = obj.Basis;      % Basis to be used
    
    % Preallocation of the control points and the expansion polynomials
    P = zeros(length(n), max(n)+1); 
    B = obj.state_basis(n, L, basis, Grid);

    % Compute the position control points leveraging the complete state vector
    if (L > 1)
        C = [s(1:Problem.StateDim,:) s(1+Problem.StateDim:2*Problem.StateDim,:)];
    else
        C = s(1:Problem.StateDim,:);       
    end

    for i = 1:length(n)
        if (L > 1)
            A = [B{i}(1:n(i)+1,:) B{i}(n(i)+2:2*(n(i)+1),:)];
        else
            A = B{i}(1:n(i)+1,:);
        end
        P(i,1:n(i)+1) = C(i,:)*pinv(A);
    end
    
    % Evaluate the state vector
    C = obj.evaluate_state(n, L, P, B);
end