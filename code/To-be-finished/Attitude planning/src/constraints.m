%% Project: Shape-based attitude planning %%
% Date: 30/01/2022

%% Constraints %% 
% Function to compute the residual vector of the constraints of the problem

% Inputs: - structure system, containing the definition of the rigid body
%           of interest
%         - scalar T, the maximum torque allowed for the spacecraft
%         - vector n, the approximation degree to each position coordinate
%         - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - vector initial, the initial boundary conditions of the
%           trajectory 
%         - vector final, the initial boundary conditions of the
%           trajectory
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - vector x, the degree of freedom to be optimized 
%         - cell array B, the polynomial basis to be used
%         - string basis, the polynomial basis to be used
%         - string dynamics, the dynamics formulation to be used

% Outputs: - inequality constraint residual vector c
%          - equality constraint residual vector ceq

function [c, ceq] = constraints(system, T, initial, final, n, x, B, basis, dynamics)
    % Extract the optimization variables
    P = reshape(x(1:end-1), [length(n), max(n)+1]);     % Control points
    tf = x(end);                                        % Maneuver time

    % Constants 
    I = system.Inertia;         % Inertia dyadic of the system 
    V = system.V;               % Preference direction 
    N = system.Prohibited;      % Prohibited attitude directions
    tol = system.Tol;           % Tolerance to the protected areas

    % Boundary conditions points
    P = boundary_conditions(tf, n, initial, final, P, B, basis);

    % Trajectory evolution
    C = evaluate_state(P, B, n);

    % Protection areas computation
    if (size(N,2) == 1)
        N = repmat(N, 1, size(C,2));
    end

    angle = zeros(1,size(C,2));         % Angle to the prohibited directions
    for i = 1:size(C,2)
        aux = quaternion_product([0; V], quaternion_inverse(C(1:4,i)));
        aux = quaternion_product(C(1:4,i), aux);
        angle(i) = tol-acos(dot(aux(2:4), N(:,i)));
    end

    % Control input 
    u = acceleration_control(I, C, tf, dynamics);

    % Equalities 
    ceq = [];

    % Inequality (control authority)
    c = [sqrt(u(1,:).^2+u(2,:).^2+u(3,:).^2)-(tf^2*T*ones(1,size(u,2))) angle];
end