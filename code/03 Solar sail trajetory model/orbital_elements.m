function [a, i] = orbital_elements(pos, vel)

global mu

% calculates the orbital elements
% in particular for this case, the semi-major axis a 
% and the orbital inclination i
% with respect to the earth's orbital plane

% the input values of position and velocity are 
% in cilindrical coordinates 

% extract data
rho = pos(1);
theta = pos(2);
z = pos(3);
drho = vel(1);
dtheta = vel(2);
dz = vel(3);

% find cartesian position
x = rho*cos(theta);
y = rho*sin(theta);
r = [x, y, z];

% find cartesian velocity
dx = drho*cos(theta) - rho*dtheta*sin(theta);
dy = drho*sin(theta) + rho*dtheta*cos(theta);
v = [dx, dy, dz];

% calculate angular momentum vector
h = cross(r,v);

% calculate eccentricity vector
e = cross(v,h)./mu - r./norm(r);

% calculate inclination
i = acos(h(3)/norm(h));

% calculate semi-major axis
a = norm(h)^2/mu/(1-norm(e)^2);

end