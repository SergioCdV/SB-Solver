function [c, ceq] = constraints(P, P0, B, amax, mu, r0, tfapp, n)

% for fmincon: non-linear constraints
% c(x)   <= 0     inequality constraints
% ceq(x) == 0     equality constraints


% non-linear constraint on acceleration:
amag = acceleration(P, mu, B, r0, tfapp, n);
c = amag - amax;



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
