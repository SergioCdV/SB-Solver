function a_cilindrical = aceleration_cilindrical(coord_cylin_adim, alphan, sigma, kappa, ac, r0)

param = kappa*ac*r0^2/2/(rho^2+z^2);

arho = rho*cos(alphan)^2 + z*sin(alphan)*cos(alphan)*cos(sigma) + rho;
atheta = cos(alphan)*sin(alphan)*sin(sigma)*sqrt(rho^2+z^2);
az = z*cos(alphan)^2 - rho*sin(alphan)*cos(alphan)*cos(sigma) + z;

a_cilindrical = param*[arho, atheta, az];


end