% Iterative portion of De Casteljau's algorithm

function P = casteljau(points,L,t)

% Extract points to the first iteration (defined points)
P(:,:,1) = points;

for i=1:1:L
    for j=1:1:L-i
        % Calulate Bezier curve (recursive)
        P(:,j,i+1) = (1-t)*P(:,j,i) + t*P(:,j+1,i); 
    end
end