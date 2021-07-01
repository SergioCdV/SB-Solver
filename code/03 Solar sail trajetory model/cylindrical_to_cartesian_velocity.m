function [X, Y, Z] = cylindrical_to_cartesian_velocity(C)

%% convert cylindrical coordinates to cartesian
rho = C{1}(:,2);
theta = C{2}(:,2);
z = C{3}(:,2);

X = rho.*cos(theta);
Y = rho.*sin(theta);

i = 0;
X = X*cos(i) + z*sin(i);
Y = Y;
Z = z*cos(i) - X*sin(i);


end