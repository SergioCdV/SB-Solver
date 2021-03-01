% This function uses the classical implementation for calculation of
% the Bezier curve as defined by the determined points.
% The time vector represents the number of points at which the curve B will
% be calculated

function curve = B_classic(points,tvec)

% Determine order
n = length(points)-1;

% Find number of steps (time increments)
steps = length(tvec);

% Initialize variable for n-order curve
B = zeros(2,steps); 

% Calculation of bezier curve, recursively
for i=0:n
    
    % binomial coefficient (n i)
    bin = factorial(n)/factorial(i)/factorial(n-i);
    
    % Bernstein basis polynomial
    bern = bin.*tvec.^i.*(1-tvec).^(n-i);
    
    % bezier curve as added to that of the previous step
    B = B + points(:,i+1)*bern;
end

curve = B;