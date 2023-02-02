%% Project: 
% Date: 31/01/22

%% Flight time %%
% Function to estimate the time of flight

% Inputs: - vector x, the vector of decision variables 
%         - cell array B, the basis of polynomials to be used 
%         - vector n, containing the order of approximation of each phase
%           space

% Outputs: - scalar r, the final orbit radius to be maximized

function [r] = minimum_time(x,B,n,tf)
%     r = x(end) * tf;
% 
%     P = reshape(x(1:end-m-1), [2 max(n)+1]);
%     C = evaluate_state(P,B,n);
%     C(3:4,:) = C(3:4,:) / r;
%     r = 0; 
%     for i = 1:size(C,2)-1
%         r = r +  norm(C(1:2,i+1)-C(1:2,i)) / norm(C(3:4,i));
%     end

    r = x(end);
end