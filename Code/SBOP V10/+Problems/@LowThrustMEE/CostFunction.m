%% Project: Shape-based optimization for low-thrust transfers %%
% Date: 07/02/2023

%% Cost function %% 
% Function implementation of a cost function 

function [M, L] = CostFunction(obj, params, beta, t0, tf, t, s, u)
    % Parameters 
    mu = params(1); 

    % Minimum fuel 
    M = 0; 

    w = 1+s(2,:).*cos(t)+s(3,:).*sin(t);
    res = sqrt(mu*s(1,:)).*(w./s(1,:)).^2 +sqrt(s(1,:)/mu)./w.*(s(4,:).*sin(t) + s(5,:).*cos(t)) .* u(3,:);
    a = 1./res;
    L = a;
end

%% Auxiliary code 
% function [cost] = minimum_time(params, beta, t0, tf, s, u)
%     % Compute the longitude evolution 
%     theta0 = initial(end)-beta(1)*L(1);
%     L = theta0+thetaf*L;
% 
%     % Kinematic constraint 
%     w = 1+s(2,:).*cos(L)+s(3,:).*sin(L);
%     a = sqrt(mu*s(1,:)).*(w./s(1,:)).^2;
%     for i = 1:size(C,2)
%         B = control_input(mu, s(:,i)); 
%         a(i) = 1/(a(i)+B(6,3)*u(3,i));
%     end
% end