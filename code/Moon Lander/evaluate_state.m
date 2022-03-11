%% Project: 
% Date: 03/02/22

%% Evaluate state %%
% Function to compute the acceleration vector norm from cylindrical coordinates

% Inputs: - array P, the set of control points to estimate the position vector 
%         - array B, the polynomial basis in use in the approximation
%         - vector n, with the degrees of approximation of each position
%           coordinate

% Outputs: - array C, the 9x1 state vector 

function [C] = evaluate_state(P, B, n)
    % Extract the spacecraft coodinates evaluating the Bézier curve approximation
    C = zeros(3,size(B{1},2));     % Preallocation for speed

    for i = 1:size(P,1)
        % State vector fitting
        k = n(i)+1;
        for j = 1:3
            C(i+1*(j-1),:) = P(i,1:n(i)+1)*B{i}(1+k*(j-1):k*j,:);
        end
    end
end