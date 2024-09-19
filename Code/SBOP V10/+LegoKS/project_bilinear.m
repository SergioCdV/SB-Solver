

function [x] = project_bilinear(u, X)
    x = X; 
    r = dot(u, u, 1);

    for i = 1:size(u,2)
        L = LegoKS.KSmatrix(u(:,i));
        gamma = dot(L(:,1), X(:,i));
        alpha = dot(L(:,2), X(:,i));
        beta =  dot(L(:,3), X(:,i));
        x(:,i) = alpha * L(:,2) + beta * L(:,3) + gamma * L(:,1);
    end

    x = x ./ r;
end
  