function [P_lb, P_ub] = point_bounds(P0, n, scale_factor)

% margin is the factor by which the points are scaled

% definition of lower and upper bounds 
% P_lb = m.*ones(size(P0))/margin;
% P_ub = 2*P0(:,n+1).*ones(size(P0))*margin;

% definition of lower and upper bounds 
P_lb = min(P0')'.*ones(size(P0))/scale_factor;
P_ub = max(P0')'.*ones(size(P0))*scale_factor;

% restriction of initial and final points
P_lb(:,1) = P0(:,1);
P_ub(:,1) = P0(:,1);
P_lb(:,n+1) = P0(:,n+1);
P_ub(:,n+1) = P0(:,n+1);

end