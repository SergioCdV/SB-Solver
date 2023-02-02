function [c, ceq] = constraints(P, P0, B, amax, mu, r0, tfapp, n)

% for fmincon: non-linear constraints
% c(x)   <= 0     inequality constraints
% ceq(x) == 0     equality constraints


% non-linear constraint on acceleration:
c1 = zeros(1,m);
for i=1:m
    a_cilindrical = aceleration_cilindrical(coord_cylin_adim, alphan, sigma, kappa, ac, r0);
    aO = acceleration_components(a_cilindrical, coord_cylin_adim);
    c1(i) = acceleration_constraints(r(i), r0, aO, ac);
end

 
% % % constraints on control variables:
% % % kappa in [0,1]
% % % alphan in [0,pi/2]
% % c2 = [-kappa, kappa-1, -alphan, alphan-pi/2];
% % 
% % 
% % c = [c1, c2];

% boundary conditions: use the initial and final position to limit points
% disp(max(margin));

ceq = zeros(3,2);

for k=1:3
    L = n(k)+1;
    P0aux = P0(k,1:L);
    Paux = P(k,1:L);
    deltaP = Paux - P0aux;

    ceq(k,:) = [deltaP(1),  deltaP(L)];

end



