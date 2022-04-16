function [X, Y, Z] = cylindrical_to_cartesian(C)

%% convert cylindrical coordinates to cartesian
rho = C{1}(:,1);
theta = C{2}(:,1);
z = C{3}(:,1);

X = rho.*cos(theta);
Y = rho.*sin(theta);

i = 0;
X = X*cos(i) + z*sin(i);
Y = Y;
Z = z*cos(i) - X*sin(i);


end