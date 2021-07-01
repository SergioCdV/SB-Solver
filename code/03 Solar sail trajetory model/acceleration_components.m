function aO = acceleration_components(a_cilindrical, coord_adim)

% function to calculate the propulsive acceleration components

arho = a_cilindrical(1);
atheta = a_cilindrical(2);
az = a_cilindrical(3);

rho = coord_adim(1);
theta = coord_adim(2);
z = coord_adim(3);


axO = (arho*z-az*rho)/sqrt(rho^2+z^2);
ayO = atheta;
azO = (arho*rho+az*z)/sqrt(rho^2+z^2);

aO = [axO, ayO, azO];



end