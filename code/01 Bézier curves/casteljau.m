%% Castlejau
% Iterative portion of De Casteljau's algorithm

% Inputs: - array points, the Bézier curve control points 
%         - scalar L, the number of considered control points
%         - scalar t, the control parameter t vector 

% Output: - the resulting Bézier curve point P

function [P] = casteljau(points,L,t)
    
    % Extract points to the first iteration (defined points)
    P(:,:,1) = points;
    
    for i=1:1:L
        for j=1:1:L-i
            % Calulate Bezier curve (recursive)
            P(:,j,i+1) = (1-t)*P(:,j,i) + t*P(:,j+1,i); 
        end
    end
end