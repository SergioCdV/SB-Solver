%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 02/08/22

%% Sundman transformation 
% Compute the derivative of time with respect to the generalized anomaly

% Inputs: - string basis, the polynomial basis to be used
%         - vector n, the vector of degrees of approximation of the state
%           variables
%         - array P, the control points of the trajectory 
%         - scalar t, the value at which the trajectory shall be evaluated

% Outputs: - scalar dt, the Sundman transformation value

function [dt] = Sundman_transformation(basis, n, P, t, ~)
    B = state_basis(n,t,basis);
    C = evaluate_state(P,B,n);
    dt = sqrt(C(1,:).^2+C(3,:).^2);
end